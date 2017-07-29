#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <boost/filesystem/operations.hpp>
#include <chrono>
#include <boost/filesystem/path.hpp>
#include "EntapInit.h"
#include "ExceptionHandler.h"
#include "pstream.h"
#include "boost/program_options.hpp"
#include "EntapConsts.h"
#include "EntapExecute.h"

namespace boostPO = boost::program_options;
namespace Chrono = std::chrono;

enum States {
    PARSE_ARGS            = 0x01,
    INIT_ENTAP            = 0x02,
    INIT_ENTAP_SUCCESS    = 0x04,
    EXECUTE_ENTAP         = 0x08,
    EXECUTE_ENTAP_SUCCESS = 0x16
};

bool check_key(std::string&);
std::unordered_map<std::string,std::string> parse_config(std::string&);
void print_msg(std::string);
void print_user_input(boostPO::variables_map&);
void init_log();
boostPO::variables_map parse_arguments_boost(int,const char**);
void state_summary(States);
std::string get_exe_path();
void generate_config(std::string&);

States state;   // init
std::string _outpath,_exe_path;
Chrono::time_point<Chrono::system_clock> _start_time, _end_time;

int main(int argc, const char** argv) {
    init_log();
    // TODO fix, not portable
    _exe_path = get_exe_path();
    try {
        state = PARSE_ARGS;
        std::unordered_map<std::string,std::string> config_map;
        boostPO::variables_map inputs = parse_arguments_boost(argc,argv);
        boost::filesystem::path working_dir(boost::filesystem::current_path());
        _outpath = working_dir.string() + "/" + inputs["tag"].as<std::string>() + "/";
        print_user_input(inputs);
        config_map = parse_config(_exe_path);
        if (state == INIT_ENTAP) {
            entapInit::init_entap(inputs, _exe_path, config_map);
            state = INIT_ENTAP_SUCCESS;
        } else if (state == EXECUTE_ENTAP) {
            entapExecute::execute_main(inputs, _exe_path,config_map);
            state = EXECUTE_ENTAP_SUCCESS;
        } else {
            print_msg("Error in parsing input data");
            return 1;
        }
    } catch (ExceptionHandler &e) {
        if (e.getErr_code()==ENTAP_ERR::E_SUCCESS) return 0;
        e.print_msg();
        state_summary(state);
        return 1;
    }
    _end_time = Chrono::system_clock::now();
    return 0;
}

boostPO::variables_map parse_arguments_boost(int argc, const char** argv) {
    std::string err_msg;
    std::unordered_map<std::string, std::string> input_map;
    print_msg("Parsing user input...");
    std::string input_file, exe_state, align_path, species;
    std::vector<std::string> contam_vec, ncbi_data, uniprot_data, data_path;
    float fpkm;
    // TODO do not change .entp filename warning
    try {
        boostPO::options_description description("Options");
        // TODO separate out into main options and additional config file with defaults
        description.add_options()
            ("help,h",
                 "Print help options")
            ("config",
                 "Configure enTAP for execution later (complete this step first)")
            (ENTAP_CONFIG::INPUT_FLAG_RUNPROTEIN.c_str(),
                 "Execute enTAP functionality with input protein sequences\n"
                 "(this option will skip Frame Selection portion of pipeline)")
            (ENTAP_CONFIG::INPUT_FLAG_RUNNUCLEOTIDE.c_str(),
                 "Execute enTAP functionality with input nucleotide sequences")
            ("ncbi,N",
                 boostPO::value<std::vector<std::string>>(&ncbi_data)->multitoken()
                 ->default_value(std::vector<std::string>{ENTAP_CONFIG::INPUT_UNIPROT_NULL},""),
                 "Select which NCBI database you would like to download"
                 "\nref - RefSeq database")
            (ENTAP_CONFIG::INPUT_FLAG_INTERPRO.c_str(),
             boostPO::value<std::vector<std::string>>()->multitoken()
                     ->default_value(std::vector<std::string>{ENTAP_CONFIG::INTERPRO_DEFAULT},""),
             "Select which protein databases you would like to download if using Interpro")
            ("uniprot,U",
                 boostPO::value<std::vector<std::string>>(&uniprot_data)->multitoken()
                 ->default_value(std::vector<std::string>{ENTAP_CONFIG::INPUT_UNIPROT_NULL},""),
                 "Select which Uniprot database you would like to download"
                 "\n100 - UniRef100...")
            (ENTAP_CONFIG::INPUT_FLAG_ONTOLOGY.c_str(),
                 boostPO::value<short>()->default_value(ENTAP_EXECUTE::EGGNOG_INT_FLAG),
                 "Specify ontology software to use\n0 - eggnog\n1 - interproscan")
            ("tag",
                 boostPO::value<std::string>()->default_value(ENTAP_EXECUTE::OUTFILE_DEFAULT),
                 "Specify species or unique tag you would like files to be saved as")
            ("database,d",
                 boostPO::value<std::vector<std::string>>(&data_path)->multitoken(),
                 "Provide the path to a separate database, however this "
                 "may prohibit taxonomic filtering.")
            (ENTAP_CONFIG::INPUT_FLAG_GO_LEVELS.c_str(),
                 boostPO::value<std::vector<short>>()->multitoken()
                 ->default_value(std::vector<short>{0,3,4},""),
                 "Gene ontology levels you would like outputted.")
            ("fpkm",
                 boostPO::value<float>(&fpkm)->default_value(ENTAP_EXECUTE::RSEM_FPKM_DEFAULT),
                 "FPKM cutoff value")
            ("e",
                 boostPO::value<double>()->default_value(ENTAP_CONFIG::E_VALUE),"Specify an e-value")
            ("version,v",
                 "Display version number")
            ("paired-end",
                 "Flag for paired end reads")
            ("threads,t",
                 boostPO::value<int>()->default_value(1),"Number of threads")
            ("align,a",
                 boostPO::value<std::string>(&align_path),"Path to BAM/SAM file")
            ("contam,c",
                 boostPO::value<std::vector<std::string>>(&contam_vec)->multitoken(),
                 "Contaminant selection")
            (ENTAP_CONFIG::INPUT_FLAG_QCOVERAGE.c_str(),
                 boostPO::value<double>()->default_value(ENTAP_CONFIG::DEFAULT_QCOVERAGE),
                 "Select minimum query coverage to be kept for similarity searching")
            (ENTAP_CONFIG::INPUT_FLAG_TCOVERAGE.c_str(),
                 boostPO::value<double>()->default_value(ENTAP_CONFIG::DEFAULT_QCOVERAGE),
                 "Select minimum target coverage to be kept for similarity searching")
            (ENTAP_CONFIG::INPUT_FLAG_SPECIES.c_str(),
                 boostPO::value<std::string>(&species),"The type of species you are analyzing if you would like"
                 "filtering based upon this separated by a '_'.\nExample: homo_sapiens")
            ("state",
                 boostPO::value<std::string>(&exe_state)->default_value("+"),
                 "Select a state value, *EXPERIMENTAL*\n""These commands will run certain "
                 "elements of the pipeline and stop at certain locations, as such"
                 "there are several runs that may be invalid as they rely on data from another portion.\n"
                 "Examples:\n+2x Will start the pipeline from Frame selection and will run RSEM then filter the"
                 "transcriptome. It will then stop execution there specified by the x.")
            ("input,i",
                 boostPO::value<std::string>(&input_file), "Input transcriptome file")
            (ENTAP_CONFIG::INPUT_FLAG_COMPLETE.c_str(),
                 "Select this option if you have all complete proteins.\n"
                         "Note: This assumes a protein input")
            (ENTAP_CONFIG::INPUT_FLAG_OVERWRITE.c_str(),
                "Select this option if you wish to overwrite pre-existing files");
        boostPO::variables_map vm;
        //TODO verify state commands

        try {
            boostPO::store(boostPO::command_line_parser(argc,argv).options(description)
                .run(),vm);
            boostPO::notify(vm);

            if (vm.count("help")) {
                std::cout << description<<std::endl<<std::endl;
                throw(ExceptionHandler("",ENTAP_ERR::E_SUCCESS));
            }
            if (vm.count("version")) {
                std::cout<<"enTAP version: "<<ENTAP_CONFIG::ENTAP_VERSION<<std::endl;
                throw(ExceptionHandler("",ENTAP_ERR::E_SUCCESS));
            }

            bool is_config = (bool) vm.count("config");     // ignore 'config config'
            bool is_protein, is_nucleotide;
            is_protein = (bool)vm.count(ENTAP_CONFIG::INPUT_FLAG_RUNPROTEIN);
            is_nucleotide = (bool)vm.count(ENTAP_CONFIG::INPUT_FLAG_RUNNUCLEOTIDE);
            if (is_protein && is_nucleotide) {
                throw ExceptionHandler("Cannot specify both protein and nucleotide",
                    ENTAP_ERR::E_INPUT_PARSE);
            }
            bool is_run = is_protein || is_nucleotide;
            if (!is_config && !is_run) {
                err_msg = "Either config option or run option are required";
                throw(ExceptionHandler(err_msg.c_str(),ENTAP_ERR::E_INPUT_PARSE));
            }
            if (is_config && is_run) {
                throw(ExceptionHandler("Cannot specify both config and run flags",
                                       ENTAP_ERR::E_INPUT_PARSE));
            }

            if (!bool(vm.count("input")) && is_run) {
                throw(ExceptionHandler("Must enter a valid transcriptome",ENTAP_ERR::E_INPUT_PARSE));
            }

            if (ncbi_data.size() + uniprot_data.size() + data_path.size() > 5) {
                // TODO fix for certain cases like -N -N -d null
                throw ExceptionHandler("Too many databases selected, 3 is the max",
                                       ENTAP_ERR::E_INPUT_PARSE);
            }
            bool ncbi_check = true;
            for (auto const& entry: ncbi_data) {
                if (entry.compare(ENTAP_CONFIG::NCBI_REFSEQ_COMP)==0)continue;
                if (entry.compare(ENTAP_CONFIG::NCBI_NONREDUNDANT)==0)continue;
                if (entry.compare(ENTAP_CONFIG::NCBI_NULL)==0)continue;
                if (entry.compare(ENTAP_CONFIG::NCBI_REFSEQ_PLANT)==0)continue;
                ncbi_check = false;
            }
            if (!ncbi_check) {
                throw ExceptionHandler("Not a valid NCBI database",ENTAP_ERR::E_INPUT_PARSE);
            }
            bool uniprot_check = true;
            for (auto const& entry: uniprot_data) {
                if (entry.compare(ENTAP_CONFIG::INPUT_UNIPROT_SWISS)==0)continue;
                if (entry.compare(ENTAP_CONFIG::INPUT_UNIPROT_TREMBL)==0)continue;
                if (entry.compare(ENTAP_CONFIG::INPUT_UNIPROT_UR90)==0)continue;
                if (entry.compare(ENTAP_CONFIG::NCBI_NULL)==0)continue;
                uniprot_check = false;
                }
            if (!uniprot_check) {
                throw ExceptionHandler("Not a valid Uniprot database",ENTAP_ERR::E_INPUT_PARSE);
            }
            if (is_run && !vm.count("input")) {
                throw ExceptionHandler("Missing input transcriptome file", ENTAP_ERR::E_INPUT_PARSE);
            }
            if (!species.empty()) {
                if (species.find("_") == std::string::npos) {
                    throw ExceptionHandler("Invalid format of species, must be "
                         "separated by a '_'", ENTAP_ERR::E_INPUT_PARSE);
                }
            }
            if (is_config) {
                state = INIT_ENTAP;
            } else state = EXECUTE_ENTAP;

            print_msg("Success!");
            // TODO parse state option
            return vm;
        } catch (boost::program_options::required_option& e) {
            std::cout<<"Required Option"<<std::endl;
        }
    }catch (boost::program_options::error& e){
        // Unknown input
        throw ExceptionHandler(e.what(),ENTAP_ERR::E_INPUT_PARSE);
    }
}

std::unordered_map<std::string,std::string> parse_config(std::string &exe) {
    print_msg("Parsing configuration file...");
    std::unordered_map<std::string,std::string> config_map;
    std::string config_path = exe + "/" + ENTAP_CONFIG::CONFIG_FILE;
    if (!entapInit::file_exists(config_path)){
        print_msg("Config file not found, generating new file...");
        try {
            generate_config(config_path);
        } catch (std::exception &e){
                throw ExceptionHandler(e.what(),ENTAP_ERR::E_CONFIG_CREATE);
        }
        print_msg("Config file successfully created");
    }
    print_msg("Config file found at: " + config_path);
    std::ifstream in_file(config_path);
    std::string line,key;
    while (std::getline(in_file,line)) {
        std::istringstream in_line(line);
        if (std::getline(in_line,key,'=')) {
            if (!check_key(key)) {
                throw ExceptionHandler("Incorrect format in config file",
                ENTAP_ERR::E_CONFIG_PARSE);
            }
            std::string val;
            if (std::getline(in_line,val)) {
                if (val.size()<=1) val = "";
                config_map.emplace(key,val);
            }
        }
    }
    print_msg("Success!");
    return config_map;
}

void generate_config(std::string &path) {
    std::ofstream config_file(ENTAP_CONFIG::CONFIG_FILE, std::ios::out | std::ios::app);
    config_file <<
                ENTAP_CONFIG::KEY_UNIPROT_SWISS             +"=\n"+
                ENTAP_CONFIG::KEY_UNIPROT_TREMBL            +"=\n"+
                ENTAP_CONFIG::KEY_UNIPROT_UR90              +"=\n"+
                ENTAP_CONFIG::KEY_UNIPROT_UR100             +"=\n"+
                ENTAP_CONFIG::KEY_NCBI_NR                   +"=\n"+
                ENTAP_CONFIG::KEY_NCBI_REFSEQ_COMPLETE      +"=\n"+
                ENTAP_CONFIG::KEY_NCBI_REFSEQ_SEPARATE      +"=\n"+
                ENTAP_CONFIG::KEY_DIAMOND_EXE               +"=\n"+
                ENTAP_CONFIG::KEY_RSEM_EXE                  +"=\n"+
                ENTAP_CONFIG::KEY_GENEMARK_EXE              +"=\n"+
                ENTAP_CONFIG::KEY_EGGNOG_EXE                +"=\n"+
                ENTAP_CONFIG::KEY_INTERPRO_EXE
                << std::endl;
    config_file.close();
}

bool check_key(std::string& key) {
    if (key.compare(ENTAP_CONFIG::KEY_NCBI_NR)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_NCBI_REFSEQ_COMPLETE)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_NCBI_REFSEQ_SEPARATE)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_UNIPROT_SWISS)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_UNIPROT_TREMBL)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_UNIPROT_UR100)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_DIAMOND_EXE)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_GENEMARK_EXE)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_UNIPROT_UR90)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_EGGNOG_EXE)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_INTERPRO_EXE)==0) return true;
    return key.compare(ENTAP_CONFIG::KEY_RSEM_EXE) == 0;
}

void init_log() {
    remove("debug.txt");
    entapInit::print_msg("Start - enTAP");
}

void print_msg(std::string msg) {
    Chrono::time_point<Chrono::system_clock> current = Chrono::system_clock::now();
    std::time_t time = Chrono::system_clock::to_time_t(current);
    std::string out_time(std::ctime(&time));
    std::ofstream log_file("debug.txt", std::ios::out | std::ios::app);
    log_file << out_time.substr(0,out_time.length()-1) << ": " + msg << std::endl;
    log_file.close();
}

void print_user_input(boostPO::variables_map &map) {
    remove(std::string(_outpath + ENTAP_CONFIG::LOG_FILENAME).c_str());
    boost::filesystem::create_directories(_outpath);
    std::stringstream ss;
    ss << ENTAP_STATS::SOFTWARE_BREAK <<
          "enTAP Run Information\n"   <<
          ENTAP_STATS::SOFTWARE_BREAK;
    _start_time = Chrono::system_clock::now();
    std::time_t time = Chrono::system_clock::to_time_t(_start_time);
    ss << "Current enTAP Version: " << ENTAP_CONFIG::ENTAP_VERSION  <<
          "\nStart time: "          << std::ctime(&time)            <<
          "\nYour working directory has been set to: "  << _outpath <<
          "\nYour execution directory has been set to: "<< _exe_path<<'\n';

    for (const auto& it : map) {
        std::string key = it.first.c_str();
        ss << "\n" << key << ": ";
        auto& value = it.second.value();
        if (auto v = boost::any_cast<std::string>(&value)) {
            ss << *v;
        } else if (auto v = boost::any_cast<std::vector<std::string>>(&value)) {
            if (v->size()>0) {
                for (auto const& val:*v) {
                    ss << val << " ";
                }
            } else ss << "null";
        } else if (auto v = boost::any_cast<float>(&value)){
            ss << *v;
        } else if (auto v = boost::any_cast<double>(&value)) {
            ss << *v;
        } else if (auto v = boost::any_cast<int>(&value)) {
            ss << *v;
        } else if (auto v = boost::any_cast<std::vector<short>>(&value)) {
            for (auto const& val:*v) {
                ss << val << " ";
            }
        } else ss << "null";
    }
    std::string output = ss.str() + "\n";
    entapExecute::print_statistics(output,_outpath);
    print_msg(output+"\n");
}

std::string get_exe_path() {
    //TODO check different systems
    char buff[1024];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
        buff[len] = '\0';
        std::string path = std::string(buff);
        boost::filesystem::path p(path);p.remove_filename();
        return p.string();
    }
    return "";
}

void state_summary(States st) {
    switch (st) {
        case(PARSE_ARGS):
            break;
        default:
            break;
    }
}