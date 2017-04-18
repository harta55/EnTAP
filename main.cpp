#include <iostream>
#include <fstream>
#include <array>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <boost/filesystem/operations.hpp>
#include "EntapInit.h"
#include "ExceptionHandler.h"
#include "pstream.h"
#include "boost/program_options.hpp"
#include "EntapConsts.h"
#include "EntapExecute.h"

namespace boostPO = boost::program_options;

enum States {
    PARSE_ARGS          = 0x01,
    INIT_ENTAP          = 0x02,
    INIT_ENTAP_SUCCESS  = 0x04,
    EXECUTE_ENTAP       = 0x08
};

void print_msg(std::string);
void init_log();
boostPO::variables_map parse_arguments_boost(int,const char**);
void state_summary(States);
std::string get_exe_path();

States state;   // init

int main(int argc, const char** argv) {
    init_log();
    // TODO fix, not portable
    std::string exe_path = get_exe_path();
    try {
        state = PARSE_ARGS;
        boostPO::variables_map inputs = parse_arguments_boost(argc,argv);
        if (state == INIT_ENTAP) {
            entapInit::init_entap(inputs, exe_path);  // todo state input 1x, user wants to start at 1 and stop
        } else if (state == EXECUTE_ENTAP) {
            entapExecute::execute_main(inputs, exe_path);
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
    return 0;
}

boostPO::variables_map parse_arguments_boost(int argc, const char** argv) {
    std::string err_msg;
    std::unordered_map<std::string, std::string> input_map;
    print_msg("Parsing user input...");
    std::string input_file, exe_state, align_path;
    std::vector<std::string> contam_vec, ncbi_data, uniprot_data, data_path;
    float fpkm;
    // TODO do not change .entp filename warning
    // TODO specify an output title (species) to name everything
    try {
        boostPO::options_description description("Options");
        // TODO separate out into main options and additional config file with defaults
        description.add_options()
                ("help,h", "Print help options")
                ("config", "Configure enTAP for execution later (complete this step first)")
                ("run", "Execute enTAP functionality")
                ("ncbi,N", boostPO::value<std::vector<std::string>>(&ncbi_data)->multitoken()
                        ->default_value(std::vector<std::string>{ENTAP_CONFIG::NCBI_DEFAULT},""),"Select which NCBI database you would like to download"
                        "\nref - RefSeq database")
                ("uniprot,U", boostPO::value<std::vector<std::string>>(&uniprot_data)->multitoken()
                         ->default_value(std::vector<std::string>{ENTAP_CONFIG::INPUT_UNIPROT_DEFAULT},""),
                        "Select which Uniprot database you would like to download"
                        "\n100 - UniRef100...")
                //multiple entries
                ("database,d", boostPO::value<std::vector<std::string>>(&data_path)->multitoken(),
                        "Provide the path to a separate database, however this "
                        "may prohibit taxonomic filtering.")
                ("fpkm",boostPO::value<float>(&fpkm)->default_value(ENTAP_EXECUTE::RSEM_FPKM_DEFAULT),
                 "FPKM cutoff value")
                ("e",boostPO::value<double>()->default_value(ENTAP_CONFIG::E_VALUE),"Specify an e-value")
                ("version,v", "Display version number")
                ("paired-end","Flag for paired end reads")
                ("threads,t",boostPO::value<int>()->default_value(1),"Number of threads")
                ("align,a", boostPO::value<std::string>(&align_path),"Path to BAM/SAM file")
                ("contam,c", boostPO::value<std::vector<std::string>>(&contam_vec)->multitoken(),"Contaminant selection")
                ("state", boostPO::value<std::string>(&exe_state)->default_value("+"),"Select a state value")
                ("input,i",boostPO::value<std::string>(&input_file), "Input transcriptome file");
        boostPO::variables_map vm;

        try {
            boostPO::store(boostPO::command_line_parser(argc,argv).options(description)
                .run(),vm);
            boostPO::notify(vm);

            if (vm.count("help")) {
                std::cout << description<<std::endl<<std::endl;
                throw(ExceptionHandler("",ENTAP_ERR::E_SUCCESS));
            }
            if (vm.count("version")) {
                std::cout<<"enTAP version 0.1.0"<<std::endl;
                throw(ExceptionHandler("",ENTAP_ERR::E_SUCCESS));
            }
            if (!bool(vm.count("input"))) {
                throw(ExceptionHandler("Must enter a valid transcriptome",ENTAP_ERR::E_INPUT_PARSE));
            }
            bool is_config = (bool) vm.count("config");     // ignore 'config config'
            bool is_run = (bool) vm.count("run");

            if (!is_config && !is_run) {
                err_msg = "Either config option or run option are required";
                throw(ExceptionHandler(err_msg.c_str(),ENTAP_ERR::E_INPUT_PARSE));
            }
            if (is_config && is_run) {
                throw(ExceptionHandler("Cannot specify both config and run flags",
                                       ENTAP_ERR::E_INPUT_PARSE));
            }

            if (ncbi_data.size() + uniprot_data.size() + data_path.size() > 5) {
                // TODO fix for certain cases like -N -N -d null
                throw ExceptionHandler("Too many databases selected, 3 is the max",
                                       ENTAP_ERR::E_INPUT_PARSE);
            }
            bool ncbi_check = false;
            for (auto const& entry: ncbi_data) {
                if (entry.compare(ENTAP_CONFIG::NCBI_REFSEQ_COMP)==0){
                    ncbi_check=true; break;
                }
                if (entry.compare(ENTAP_CONFIG::NCBI_NONREDUNDANT)==0) {
                    ncbi_check=true;break;
                }
                if (entry.compare(ENTAP_CONFIG::NCBI_NULL)==0){
                    ncbi_check=true;break;
                }
                if (entry.compare(ENTAP_CONFIG::NCBI_REFSEQ_PLANT)==0){
                    ncbi_check=true;break;
                }
            }
            if (!ncbi_check) {
                throw ExceptionHandler("Not a valid NCBI database",ENTAP_ERR::E_INPUT_PARSE);
            }
            bool uniprot_check = false;
            for (auto const& entry: uniprot_data) {

                if (entry.compare(ENTAP_CONFIG::INPUT_UNIPROT_SWISS)==0){
                    uniprot_check=true; break;
                }
                if (entry.compare(ENTAP_CONFIG::INPUT_UNIPROT_TREMBL)==0) {
                    uniprot_check=true;break;
                }
                if (entry.compare(ENTAP_CONFIG::INPUT_UNIPROT_UR90)==0){
                    uniprot_check=true;break;
                }
                if (entry.compare(ENTAP_CONFIG::NCBI_NULL)==0){
                    uniprot_check=true;break;
                }
            }
            if (!uniprot_check) {
                throw ExceptionHandler("Not a valid Uniprot database",ENTAP_ERR::E_INPUT_PARSE);
            }
            // TODO check unknown database

            if (is_run && !vm.count("input")) {
                throw ExceptionHandler("Missing input transcriptome file", ENTAP_ERR::E_INPUT_PARSE);
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

void init_log() {
    remove("debug.txt");
    print_msg("Start - enTAP");
}

void print_msg(std::string msg) {
    time_t rawtime;
    time(&rawtime);
    std::string date_time = ctime(&rawtime);
    std::ofstream log_file("debug.txt", std::ios::out | std::ios::app);
    log_file << date_time.substr(0, date_time.size() - 2)
                 + ": " + msg << std::endl;
    log_file.close();
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
    return "";      // some error going on with stream
    boost::filesystem::path p(buff);p.parent_path();
    std::cout<<p.string()<<std::endl;
}

void state_summary(States st) {
    switch (st) {
        case(PARSE_ARGS):
            break;
        default:
            break;
    }
}