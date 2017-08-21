//
// Created by harta on 8/4/17.
//


//*********************** Includes *****************************
#include <unordered_map>
#include <boost/program_options/options_description.hpp>
#include <iostream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <chrono>
#include "UserInput.h"
#include "EntapGlobals.h"
#include "ExceptionHandler.h"

//**************************************************************

bool _config;


/**
 * ======================================================================
 * Function std::pair<bool, boostPO::variables_map> entap_user_parse
 *                              (int            argc,
 *                               const char**   argv)
 *
 * Description          - Handles overall parsing and verifying of user
 *                        flags
 *                      - Flags parsed through Boost program options
 *
 * Notes                - Called from main
 *
 * @param argc          - Pushed from main
 * @param argv          - Pushed from main
 * @return              - Pair containing flag if configuration is selected
 *                        and map of user parsed flags
 * ======================================================================
 */
std::pair<bool, boostPO::variables_map> entap_user_parse(int argc, const char** argv) {

    boostPO::variables_map user_inputs;
    std::pair<bool, boostPO::variables_map> out_pair;

    try {
        _config = false;
        user_inputs = parse_arguments_boost(argc, argv);
        verify_user_input(user_inputs);

    } catch (const ExceptionHandler &e) {throw e;}
    out_pair.first  = _config;
    out_pair.second = user_inputs;
    return out_pair;
};


/**
 * ======================================================================
 * Function boostPO::variables_map parse_arguments_boost
 *                              (int            argc,
 *                               const char**   argv)
 *
 * Description          - Utilizes boost libraries to parse user input
 *                        arguments
 *
 * Notes                - None
 *
 * @param argc          - Pushed from main
 * @param argv          - Pushed from main
 * @return              - Variable map of user input flags
 * ======================================================================
 */
boost::program_options::variables_map parse_arguments_boost(int argc, const char** argv) {
    print_debug("Parsing user input...");

    std::unordered_map<std::string, std::string> input_map;

    try {
        boostPO::options_description description("Options");
        // TODO separate out into main options and additional config file with defaults
        description.add_options()
                ("help,h",
                 "Print help options")
                (ENTAP_CONFIG::INPUT_FLAG_CONFIG.c_str(),
                 "Configure EnTAP for execution later (complete this step first)")
                (ENTAP_CONFIG::INPUT_FLAG_RUNPROTEIN.c_str(),
                 "Execute EnTAP functionality with input protein sequences\n"
                         "(this option will skip Frame Selection portion of pipeline)")
                (ENTAP_CONFIG::INPUT_FLAG_RUNNUCLEOTIDE.c_str(),
                 "Execute EnTAP functionality with input nucleotide sequences")
                ("ncbi,N",
                 boostPO::value<std::vector<std::string>>()->multitoken()
                         ->default_value(std::vector<std::string>{ENTAP_CONFIG::INPUT_UNIPROT_NULL},""),
                 "Coming soon!")
                (ENTAP_CONFIG::INPUT_FLAG_INTERPRO.c_str(),
                 boostPO::value<std::vector<std::string>>()->multitoken()
                         ->default_value(std::vector<std::string>{INTERPRO_DEFAULT},""),
                 "Select which protein databases you would like to download if using Interpro")
                ("uniprot,U",
                 boostPO::value<std::vector<std::string>>()->multitoken()
                         ->default_value(std::vector<std::string>{ENTAP_CONFIG::INPUT_UNIPROT_NULL},""),
                 "Coming Soon!")
                (ENTAP_CONFIG::INPUT_FLAG_ONTOLOGY.c_str(),
                 boostPO::value<short>()->default_value(ENTAP_EXECUTE::EGGNOG_INT_FLAG),
                 "Specify ontology software to use\n0 - eggnog\n1 - interproscan")
                ("tag",
                 boostPO::value<std::string>()->default_value(OUTFILE_DEFAULT),
                 "Specify species or unique tag you would like files to be saved as")
                ("database,d",
                 boostPO::value<std::vector<std::string>>()->multitoken(),
                 "Provide the path to a separate database, however this "
                         "may prohibit taxonomic filtering.")
                (ENTAP_CONFIG::INPUT_FLAG_GO_LEVELS.c_str(),
                 boostPO::value<std::vector<short>>()->multitoken()
                         ->default_value(std::vector<short>{0,3,4},""),
                 "Gene ontology levels you would like outputted.")
                (ENTAP_CONFIG::INPUT_FLAG_FPKM.c_str(),
                 boostPO::value<float>()->default_value(RSEM_FPKM_DEFAULT),
                 "FPKM cutoff value")
                ("e",
                 boostPO::value<float>()->default_value(E_VALUE),"Specify an e-value")
                ("version,v",
                 "Display version number")
                ("paired-end",
                 "Flag for paired end reads")
                ("threads,t",
                 boostPO::value<int>()->default_value(1),"Number of threads")
                ("align,a",
                 boostPO::value<std::string>(),"Path to BAM/SAM file")
                ("contam,c",
                 boostPO::value<std::vector<std::string>>()->multitoken(),
                 "Contaminant selection")
                (ENTAP_CONFIG::INPUT_FLAG_QCOVERAGE.c_str(),
                 boostPO::value<float>()->default_value(DEFAULT_QCOVERAGE),
                 "Select minimum query coverage to be kept for similarity searching")
                (ENTAP_CONFIG::INPUT_FLAG_EXE_PATH.c_str(),
                 boostPO::value<std::string>(),
                 "Specify path to EnTAP exe if it is not detected by the program.")
                (ENTAP_CONFIG::INPUT_FLAG_DATA_OUT.c_str(),
                 boostPO::value<std::string>(),
                 "Specify output directory for diamond formatted databases.")
                (ENTAP_CONFIG::INPUT_FLAG_TCOVERAGE.c_str(),
                 boostPO::value<float>()->default_value(DEFAULT_TCOVERAGE),
                 "Select minimum target coverage to be kept for similarity searching")
                (ENTAP_CONFIG::INPUT_FLAG_SPECIES.c_str(),
                 boostPO::value<std::string>(),"The type of species you are analyzing if you would like"
                         "filtering based upon this separated by a '_'.\nExample: homo_sapiens")
                ("state",
                 boostPO::value<std::string>()->default_value("+"),
                 "Select a state value, *EXPERIMENTAL*\n""These commands will run certain "
                         "elements of the pipeline and stop at certain locations, as such"
                         "there are several runs that may be invalid as they rely on data from another portion.\n"
                         "Examples:\n+2x Will start the pipeline from Frame selection and will run RSEM then filter the"
                         "transcriptome. It will then stop execution there specified by the x.")
                ("input,i",
                 boostPO::value<std::string>(), "Input transcriptome file")
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
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_VERSION)) {
                std::cout<<"enTAP version: "<<ENTAP_CONFIG::ENTAP_VERSION<<std::endl;
                throw(ExceptionHandler("",ENTAP_ERR::E_SUCCESS));
            }
            print_debug("Success!");
            return vm;
        } catch (boost::program_options::required_option& e) {
            throw ExceptionHandler(e.what(),ENTAP_ERR::E_INPUT_PARSE);
        }
    }catch (boost::program_options::error& e){
        // Unknown input
        throw ExceptionHandler(e.what(),ENTAP_ERR::E_INPUT_PARSE);
    }
}


/**
 * ======================================================================
 * Function void verify_user_input(boostPO::variables_map&       vm)
 *
 * Description          - Mangages ensuring user input is valid and
 *                        will not cause issues downstream
 *
 * Notes                - None
 *
 * @param vm            - User variable map of flags
 * @return              - None
 * =====================================================================
 */
void verify_user_input(boostPO::variables_map& vm) {

    bool                     is_protein;
    bool                     is_nucleotide;
    bool                     is_config;
    bool                     is_run;
    std::string              species;
    std::string              input_tran_path;

    is_config = (bool) vm.count(ENTAP_CONFIG::INPUT_FLAG_CONFIG);     // ignore 'config config'
    is_protein = (bool)vm.count(ENTAP_CONFIG::INPUT_FLAG_RUNPROTEIN);
    is_nucleotide = (bool)vm.count(ENTAP_CONFIG::INPUT_FLAG_RUNNUCLEOTIDE);
    if (is_protein && is_nucleotide) {
        throw ExceptionHandler("Cannot specify both protein and nucleotide input flags",
                               ENTAP_ERR::E_INPUT_PARSE);
    }
    is_run = is_protein || is_nucleotide;
    if (!is_config && !is_run) {
        throw(ExceptionHandler("Either config option or run option are required",
                               ENTAP_ERR::E_INPUT_PARSE));
    }
    if (is_config && is_run) {
        throw(ExceptionHandler("Cannot specify both config and run flags",
                               ENTAP_ERR::E_INPUT_PARSE));
    }
    _config = is_config;

    try {
        verify_databases(vm);

        // Handle EnTAP execution commands
        if (is_run) {

            // Verify input transcriptome
            if (!vm.count(ENTAP_CONFIG::INPUT_FLAG_TRANSCRIPTOME)) {
                throw(ExceptionHandler("Must enter a valid transcriptome",ENTAP_ERR::E_INPUT_PARSE));
            } else {
                input_tran_path = vm[ENTAP_CONFIG::INPUT_FLAG_TRANSCRIPTOME].as<std::string>();
                if (!file_exists(input_tran_path)) {
                    throw(ExceptionHandler("Transcriptome not found at: " + input_tran_path,
                                           ENTAP_ERR::E_INPUT_PARSE));
                }
            }

            // Verify species for taxonomic relevance
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_SPECIES)) {
                verify_species(vm);
            }

            // Verify path + extension for alignment file
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_ALIGN)) {
                std::string align_file = vm[ENTAP_CONFIG::INPUT_FLAG_ALIGN].as<std::string>();
                std::string align_ext = boostFS::path(align_file).extension().string();
                std::transform(align_ext.begin(), align_ext.end(), align_ext.begin(), ::tolower);
                if (align_ext.compare(".sam") != 0 && align_ext.compare(".bam") != 0) {
                    throw ExceptionHandler("Alignment file must have a .bam or .sam extension",
                                           ENTAP_ERR::E_INPUT_PARSE);
                }
                if (!file_exists(align_file)) {
                    throw ExceptionHandler("Invalid file path for BAM/SAM file, exiting...",
                                           ENTAP_ERR::E_INIT_TAX_READ);
                }
            }
        }

    }catch (const ExceptionHandler &e) {throw e;}
}


/**
 * ======================================================================
 * Function void verify_databases(boostPO::variables_map& vm)
 *
 * Description          - Ensures the user is entering valid databases
 *                        and flags
 *
 * Notes                - Not really currently used, will be updated
 *
 * @param exe           - Boost variable map of user input
 * @return              - None
 * ======================================================================
 */
void verify_databases(boostPO::variables_map& vm) {

    databases_t     uniprot_data;
    databases_t     ncbi_data;
    databases_t     other_data;
    bool            ncbi_check;
    bool            uniprot_check;

    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_UNIPROT)) {
        uniprot_data = vm[ENTAP_CONFIG::INPUT_FLAG_UNIPROT].as<databases_t>();
    }
    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_NCBI_1)) {
        ncbi_data = vm[ENTAP_CONFIG::INPUT_FLAG_NCBI_1].as<databases_t>();
    }
    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_DATABASE)) {
        other_data = vm[ENTAP_CONFIG::INPUT_FLAG_DATABASE].as<databases_t>();
    }

    if (ncbi_data.size() + uniprot_data.size() + other_data.size() > MAX_DATABASE_SIZE) {
        // TODO fix for certain cases like -N -N -d null
        throw ExceptionHandler("Too many databases selected, 3 is the max",
                               ENTAP_ERR::E_INPUT_PARSE);
    }
    ncbi_check = true;
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
    uniprot_check = true;
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
}


/**
 * ======================================================================
 * Function std::unordered_map<std::string,std::string> parse_config(std::string &exe)
 *
 * Description          - Manages parsing and verifying the EnTAP configuration
 *                        file provided by the user
 *                      - Calls verification of each key to ensure user
 *                        has not changed keys
 *                      - Generates configuration file if one is not found
 *
 * Notes                - Entry
 *
 * @param exe           - Path to EnTAP executable or main directory
 * @return              - Map of keys to user parameters
 * ======================================================================
 */
std::unordered_map<std::string,std::string> parse_config(std::string &exe) {
    print_debug("Parsing configuration file...");

    std::unordered_map<std::string,std::string> config_map;
    std::string                                 config_path;
    std::string                                 line;
    std::string                                 key;
    std::string                                 val;

    config_path = (boostFS::path(exe) / boostFS::path(ENTAP_CONFIG::CONFIG_FILE)).string();
    if (!file_exists(config_path)){
        print_debug("Config file not found, generating new file...");
        try {
            generate_config(config_path);
        } catch (std::exception &e){
            throw ExceptionHandler(e.what(),ENTAP_ERR::E_CONFIG_CREATE);
        }
        print_debug("Config file successfully created");
    }
    print_debug("Config file found at: " + config_path);
    std::ifstream in_file(config_path);
    while (std::getline(in_file,line)) {
        std::istringstream in_line(line);
        if (std::getline(in_line,key,'=')) {
            if (!check_key(key)) {
                throw ExceptionHandler("Incorrect format in config file",
                                       ENTAP_ERR::E_CONFIG_PARSE);
            }
            if (std::getline(in_line,val)) {
                if (val.size()<=1) val = "";
                config_map.emplace(key,val);
            }
        }
    }
    print_debug("Success!");
    return config_map;
}


/**
 * ======================================================================
 * Function void generate_config(std::string        &path)
 *
 * Description          - Generates configuration file if it is not found
 *                        previously
 *
 * Notes                - Called from function parse_config()
 *
 * @param path          - Path to where the file should be written
 * @return              - None
 * =====================================================================
 */
void generate_config(std::string &path) {
    std::ofstream config_file(path, std::ios::out | std::ios::app);
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
                ENTAP_CONFIG::KEY_EGGNOG_DOWN               +"=\n"+
                ENTAP_CONFIG::KEG_EGGNOG_DMND_DB            +"=\n"+
                ENTAP_CONFIG::KEY_INTERPRO_EXE
                << std::endl;
    config_file.close();
}


/**
 * ======================================================================
 * Function bool check_key(std::string&     key)
 *
 * Description          - Ensures EnTAP configuration file has valid
 *                        entries and has not been edited by user
 *
 * Notes                - Called from function parse_config()
 *
 * @param key           - Key from configuration file
 * @return              - Flag if key is valid or not
 * =====================================================================
 */
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
    if (key.compare(ENTAP_CONFIG::KEY_EGGNOG_DOWN)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEG_EGGNOG_DMND_DB)==0) return true;
    if (key.compare(ENTAP_CONFIG::KEY_INTERPRO_EXE)==0) return true;
    return key.compare(ENTAP_CONFIG::KEY_RSEM_EXE) == 0;
}


/**
 * ======================================================================
 * Function void print_user_input(boostPO::variables_map        &map)
 *
 * Description          - Handles printing of user selected flags to
 *                        EnTAP statistics/log file
 *
 * Notes                - Called from main
 *
 * @param map           - Boost parsed map of user inputs
 * @param exe           - EnTAP exe/main directory
 * @param out           - Working directory
 * @return              - None
 *
 * =====================================================================
 */
void print_user_input(boostPO::variables_map &map, std::string& exe, std::string &out) {

    std::string         output;
    std::stringstream   ss;
    std::time_t         time;
    std::chrono::time_point<std::chrono::system_clock> _start_time;


    _start_time = std::chrono::system_clock::now();
    time = std::chrono::system_clock::to_time_t(_start_time);

    ss <<
       ENTAP_STATS::SOFTWARE_BREAK <<
       "enTAP Run Information\n"   <<
       ENTAP_STATS::SOFTWARE_BREAK <<
       "Current enTAP Version: "   << ENTAP_CONFIG::ENTAP_VERSION  <<
       "\nStart time: "            << std::ctime(&time)            <<
       "\nWorking directory has been set to: "  << out        <<
       "\nExecution directory has been set to: "<< exe        <<'\n';

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
    output = ss.str() + "\n";
    print_statistics(output);
    print_debug(output+"\n");
}


/**
 * ======================================================================
 * Function void verify_species(boostPO::variables_map &map)
 *
 * Description          - Verify species/tax level input by the user
 *                      - Ensure it can be found within the tax database
 *                      - Ensure it's in the right format
 *
 * Notes                - None
 *
 * @param exe           - Boost map of user inputs
 * @return              - None
 * ======================================================================
 */
void verify_species(boostPO::variables_map &map) {

    std::string species;

    species = map[ENTAP_CONFIG::INPUT_FLAG_SPECIES].as<std::string>();

    if (species.empty()) return;

    if (species.find("_") == std::string::npos) {
        throw ExceptionHandler("Invalid format of species, must be "
                                       "separated by a '_'", ENTAP_ERR::E_INPUT_PARSE);
    }

    // TODO check it can be found within tax database
}