/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
 *
 * This file is part of EnTAP.
 *
 * EnTAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EnTAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/


//*********************** Includes *****************************
#include <unordered_map>
#include <boost/program_options/options_description.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <chrono>
#include "UserInput.h"
#include "EntapGlobals.h"
#include "ExceptionHandler.h"
#include "GraphingManager.h"
#include "SimilaritySearch.h"
#include "common.h"
#include "FileSystem.h"
#include "ontology/ModEggnog.h"

//**************************************************************

std::string RSEM_EXE_DIR;
std::string GENEMARK_EXE;
std::string DIAMOND_EXE;
std::string EGG_EMAPPER_EXE;
std::string EGG_SQL_DB_PATH;
std::string EGG_DOWNLOAD_EXE;
std::string INTERPRO_EXE;
std::string TAX_DB_PATH;
std::string TAX_DOWNLOAD_EXE;
std::string GO_DB_PATH;
std::string GRAPHING_EXE;


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
//    FS_dprint("Parsing user input...");

    std::unordered_map<std::string, std::string> input_map;

    try {
        boostPO::options_description description("Options");
        // TODO separate out into main options and additional config file with defaults
        description.add_options()
                ("help,h", DESC_HELP)
                (ENTAP_CONFIG::INPUT_FLAG_CONFIG.c_str(),DESC_CONFIG)
                (ENTAP_CONFIG::INPUT_FLAG_RUNPROTEIN.c_str(),DESC_RUN_PROTEIN)
                (ENTAP_CONFIG::INPUT_FLAG_RUNNUCLEOTIDE.c_str(),DESC_RUN_NUCLEO)
                (ENTAP_CONFIG::INPUT_FLAG_UNINFORM.c_str(), boostPO::value<std::string>(),DESC_UNINFORMATIVE)
                ("ncbi,N",
                 boostPO::value<std::vector<std::string>>()->multitoken()
                         ->default_value(std::vector<std::string>{ENTAP_CONFIG::INPUT_UNIPROT_NULL},""),
                 DESC_COMING_SOON)
                (ENTAP_CONFIG::INPUT_FLAG_INTERPRO.c_str(),
                 boostPO::value<std::vector<std::string>>()->multitoken()
                         ->default_value(std::vector<std::string>{INTERPRO_DEFAULT},""),DESC_INTER_DATA)
                ("uniprot,U",
                 boostPO::value<std::vector<std::string>>()->multitoken()
                         ->default_value(std::vector<std::string>{ENTAP_CONFIG::INPUT_UNIPROT_NULL},""),
                 DESC_COMING_SOON)
                (ENTAP_CONFIG::INPUT_FLAG_ONTOLOGY.c_str(),
                 boostPO::value<std::vector<uint16>>()->multitoken()
                 ->default_value(std::vector<uint16>{ENTAP_EXECUTE::EGGNOG_INT_FLAG},""),DESC_ONTOLOGY_FLAG)
                (ENTAP_CONFIG::INPUT_FLAG_GRAPH.c_str(),DESC_GRAPHING)
                (ENTAP_CONFIG::INPUT_FLAG_TAG.c_str(),
                 boostPO::value<std::string>()->default_value(OUTFILE_DEFAULT),DESC_OUT_FLAG)
                ("database,d",
                 boostPO::value<std::vector<std::string>>()->multitoken(),DESC_DATABASE)
                (ENTAP_CONFIG::INPUT_FLAG_GO_LEVELS.c_str(),
                 boostPO::value<std::vector<uint16>>()->multitoken()
                         ->default_value(std::vector<uint16>{0,3,4},""), DESC_ONT_LEVELS)
                (ENTAP_CONFIG::INPUT_FLAG_FPKM.c_str(),
                 boostPO::value<fp32>()->default_value(RSEM_FPKM_DEFAULT), DESC_FPKM)
                (ENTAP_CONFIG::INPUT_FLAG_E_VAL.c_str(),
                 boostPO::value<fp32>()->default_value(E_VALUE),DESC_EVAL)
                ("version,v", "Display version number")
                (ENTAP_CONFIG::INPUT_FLAG_SINGLE_END.c_str(), DESC_SINGLE_END)
                ("threads,t",
                 boostPO::value<int>()->default_value(1),DESC_THREADS)
                ("align,a", boostPO::value<std::string>(),DESC_ALIGN_FILE)
                ("contam,c",
                 boostPO::value<std::vector<std::string>>()->multitoken(),DESC_CONTAMINANT)
                (ENTAP_CONFIG::INPUT_FLAG_TRIM.c_str(), DESC_TRIM)
                (ENTAP_CONFIG::INPUT_FLAG_QCOVERAGE.c_str(),
                 boostPO::value<fp32>()->default_value(DEFAULT_QCOVERAGE), DESC_QCOVERAGE)
                (ENTAP_CONFIG::INPUT_FLAG_EXE_PATH.c_str(), boostPO::value<std::string>(), DESC_EXE_PATHS)
                (ENTAP_CONFIG::INPUT_FLAG_DATA_OUT.c_str(), boostPO::value<std::string>(), DESC_DATA_OUT)
                (ENTAP_CONFIG::INPUT_FLAG_TCOVERAGE.c_str(),
                 boostPO::value<fp32>()->default_value(DEFAULT_TCOVERAGE), DESC_TCOVERAGE)
                (ENTAP_CONFIG::INPUT_FLAG_SPECIES.c_str(), boostPO::value<std::string>(),DESC_TAXON)
                (ENTAP_CONFIG::INPUT_FLAG_STATE.c_str(),
                 boostPO::value<std::string>()->default_value(DEFAULT_STATE), DESC_STATE)
                ("input,i", boostPO::value<std::string>(), DESC_INPUT_TRAN)
                (ENTAP_CONFIG::INPUT_FLAG_COMPLETE.c_str(), DESC_COMPLET_PROT)
                (ENTAP_CONFIG::INPUT_FLAG_NOCHECK.c_str(), DESC_NOCHECK)
                (ENTAP_CONFIG::INPUT_FLAG_OVERWRITE.c_str(), DESC_OVERWRITE);
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
                std::cout<<"EnTAP version: "<<ENTAP_CONFIG::ENTAP_VERSION<<std::endl;
                throw(ExceptionHandler("",ENTAP_ERR::E_SUCCESS));
            }
//            FS_dprint("Success!");
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
bool verify_user_input(boostPO::variables_map& vm) {

    bool                     is_interpro;
    bool                     is_protein;
    bool                     is_nucleotide;
    bool                     is_config;
    bool                     is_run;
    std::string              species;
    std::string              input_tran_path;
    std::vector<uint16>      ont_flags;

    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_GRAPH)) {
        if (!FS_file_exists(GRAPHING_EXE)) {
            std::cout<<"Graphing is NOT enabled on this system! Graphing script could not "
                    "be found at: "<<GRAPHING_EXE << std::endl;
        }
        GraphingManager gmanager = GraphingManager(GRAPHING_EXE);
        if (gmanager.is_graphing_enabled()) {
            std::cout<< "Graphing is enabled on this system!" << std::endl;
            throw ExceptionHandler("",ENTAP_ERR::E_SUCCESS);
        } else {
            std::cout<<"Graphing is NOT enabled on this system!,"
                    " ensure that you have python with the Matplotlib module installed."<<std::endl;
            throw ExceptionHandler("",ENTAP_ERR::E_SUCCESS);
        }
    }


    // ------------ Config / Run Required beyond this point ---------------- //

    is_config     = (bool) vm.count(ENTAP_CONFIG::INPUT_FLAG_CONFIG);     // ignore 'config config'
    is_protein    = (bool)vm.count(ENTAP_CONFIG::INPUT_FLAG_RUNPROTEIN);
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

    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_NOCHECK)) return is_config;

    try {
        verify_databases(vm);

        // Handle EnTAP execution commands
        if (is_run) {

            // Verify input transcriptome
            if (!vm.count(ENTAP_CONFIG::INPUT_FLAG_TRANSCRIPTOME)) {
                throw(ExceptionHandler("Must enter a valid transcriptome",ENTAP_ERR::E_INPUT_PARSE));
            } else {
                input_tran_path = vm[ENTAP_CONFIG::INPUT_FLAG_TRANSCRIPTOME].as<std::string>();
                if (!FS_file_exists(input_tran_path)) {
                    throw(ExceptionHandler("Transcriptome not found at: " + input_tran_path,
                                           ENTAP_ERR::E_INPUT_PARSE));
                } else if (FS_file_empty(input_tran_path)) {
                    throw(ExceptionHandler("Transcriptome file empty: "+ input_tran_path,
                                            ENTAP_ERR::E_INPUT_PARSE));
                } else if (!FS_check_fasta(input_tran_path)) {
                    throw(ExceptionHandler("File not in fasta format or corrupt! "+ input_tran_path,
                                           ENTAP_ERR::E_INPUT_PARSE));
                }
            }

            // Verify species for taxonomic relevance
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_SPECIES)) {
                verify_species(vm, SPECIES);
            }

            // Verify contaminant
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_CONTAM)) {
                verify_species(vm, CONTAMINANT);
            }

            // Verify path + extension for alignment file
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_ALIGN)) {
                std::string align_file = vm[ENTAP_CONFIG::INPUT_FLAG_ALIGN].as<std::string>();
                std::string align_ext = boostFS::path(align_file).extension().string();
                std::transform(align_ext.begin(), align_ext.end(), align_ext.begin(), ::tolower);
                if (align_ext.compare(SAM_EXT) != 0 && align_ext.compare(BAM_EXT) != 0) {
                    throw ExceptionHandler("Alignment file must have a .bam or .sam extension",
                                           ENTAP_ERR::E_INPUT_PARSE);
                }
                if (!FS_file_exists(align_file)) {
                    throw ExceptionHandler("Invalid file path for BAM/SAM file, exiting...",
                                           ENTAP_ERR::E_INIT_TAX_READ);
                }
            }

            // Verify FPKM
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_FPKM)) {
                fp32 fpkm = vm[ENTAP_CONFIG::INPUT_FLAG_FPKM].as<fp32>();
                if (fpkm > FPKM_MAX || fpkm < FPKM_MIN) {
                    throw ExceptionHandler("FPKM is out of range, but be between " + std::to_string(FPKM_MIN) +
                                           " and " + std::to_string(FPKM_MAX), ENTAP_ERR::E_INPUT_PARSE);
                }
            }

            // Verify query coverage
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_QCOVERAGE)) {
                fp32 qcoverage = vm[ENTAP_CONFIG::INPUT_FLAG_QCOVERAGE].as<fp32>();
                if (qcoverage > COVERAGE_MAX || qcoverage < COVERAGE_MIN) {
                    throw ExceptionHandler("Query coverage is out of range, but be between " +
                                           std::to_string(COVERAGE_MIN) +
                                           " and " + std::to_string(COVERAGE_MAX), ENTAP_ERR::E_INPUT_PARSE);
                }
            }

            // Verify target coverage
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_TCOVERAGE)) {
                fp32 qcoverage = vm[ENTAP_CONFIG::INPUT_FLAG_TCOVERAGE].as<fp32>();
                if (qcoverage > COVERAGE_MAX || qcoverage < COVERAGE_MIN) {
                    throw ExceptionHandler("Target coverage is out of range, but be between " +
                                           std::to_string(COVERAGE_MIN) +
                                           " and " + std::to_string(COVERAGE_MAX), ENTAP_ERR::E_INPUT_PARSE);
                }
            }

            // Verify for default state, may need to do a temp run of executables to verify
            if (vm[ENTAP_CONFIG::INPUT_FLAG_STATE].as<std::string>() == DEFAULT_STATE) {
                if (!FS_file_exists(TAX_DB_PATH)) {
                    throw ExceptionHandler("Taxonomic database could not be found at: " + TAX_DB_PATH +
                                            " make sure to set the path in the configuration file",
                                            ENTAP_ERR::E_INPUT_PARSE);
                }
            }

            // Verify Ontology Flags
            is_interpro = false;
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_ONTOLOGY)) {
                ont_flags = vm[ENTAP_CONFIG::INPUT_FLAG_ONTOLOGY].as<std::vector<uint16>>();
                for (int i = 0; i < ont_flags.size() ; i++) {
                    if ((ont_flags[i] > ENTAP_EXECUTE::ONTOLOGY_MAX) ||
                         ont_flags[i] < ENTAP_EXECUTE::ONTOLOGY_MIN) {
                        throw ExceptionHandler("Invalid ontology flags being used",
                                               ENTAP_ERR::E_INPUT_PARSE);
                    }
                    if (ont_flags[i] == ENTAP_EXECUTE::INTERPRO_INT_FLAG && !is_interpro) is_interpro = true;
                }
            }

            // Verify InterPro databases
            if (is_interpro) {
                if (vm.count(ENTAP_CONFIG::INPUT_FLAG_INTERPRO)) {
                    std::vector<std::string> inter_databases =
                            vm[ENTAP_CONFIG::INPUT_FLAG_INTERPRO].as<std::vector<std::string>>();
                    for (std::string &database : inter_databases) {
                        if (!verify_interpro(database)) {
                            throw ExceptionHandler("InterPro database: " + database + " invalid!",
                                                   ENTAP_ERR::E_INPUT_PARSE);
                        }
                    }
                } else {
                    throw ExceptionHandler("InterPro selected, but no databases specified.", ENTAP_ERR::E_INPUT_PARSE);
                }
            }

            // Verify uninformative file list
            if (vm.count(ENTAP_CONFIG::INPUT_FLAG_UNINFORM)) {
                std::string uninform_path = vm[ENTAP_CONFIG::INPUT_FLAG_UNINFORM].as<std::string>();
                verify_uninformative(uninform_path);
            }

            // Verify paths from state
            if (vm[ENTAP_CONFIG::INPUT_FLAG_STATE].as<std::string>().compare(DEFAULT_STATE)==0) {
                std::string state = vm[ENTAP_CONFIG::INPUT_FLAG_STATE].as<std::string>();
                // onlty handling default now
                verify_state(state, is_protein, ont_flags);
            }


        } else {
            // Must be config
            ;
        }

    }catch (const ExceptionHandler &e) {throw e;}
    return is_config;
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

    databases_t     other_data;

#if NCBI_UNIPROT
    databases_t     uniprot_data;
    databases_t     ncbi_data;
    bool            ncbi_check;
    bool            uniprot_check;

    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_UNIPROT)) {
        uniprot_data = vm[ENTAP_CONFIG::INPUT_FLAG_UNIPROT].as<databases_t>();
    }
    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_NCBI_1)) {
        ncbi_data = vm[ENTAP_CONFIG::INPUT_FLAG_NCBI_1].as<databases_t>();
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
#endif

    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_DATABASE)) {
        other_data = vm[ENTAP_CONFIG::INPUT_FLAG_DATABASE].as<databases_t>();
    }
    if (other_data.size() > MAX_DATABASE_SIZE) {
        throw ExceptionHandler("Too many databases selected, the max is " +
            std::to_string(MAX_DATABASE_SIZE), ENTAP_ERR::E_INPUT_PARSE);
    }
    for (auto const& path: other_data) {
        if (!FS_file_exists(path)) throw ExceptionHandler("Database path invalid: " +
            path, ENTAP_ERR::E_INPUT_PARSE);
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
std::unordered_map<std::string,std::string> parse_config(std::string &config,std::string &exe) {
    FS_dprint("Parsing configuration file...");

    std::unordered_map<std::string,std::string> config_map;
    std::string                                 new_config;
    std::string                                 line;
    std::string                                 key;
    std::string                                 val;

    if (!FS_file_exists(config)){
        FS_dprint("Config file not found, generating new file...");
        new_config = CONFIG_FILE;
        try {
            generate_config(new_config);
        } catch (std::exception &e){
            throw ExceptionHandler(e.what(),ENTAP_ERR::E_CONFIG_CREATE);
        }
        FS_dprint("Config file successfully created");
        throw ExceptionHandler("Configuration file generated at: " + config,
            ENTAP_ERR::E_CONFIG_CREATE_SUCCESS);
    }
    FS_dprint("Config file found at: " + config);
    std::ifstream in_file(config);
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
    FS_dprint("Success!");
    init_exe_paths(config_map,exe);
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
                KEY_DIAMOND_EXE               <<"=\n"<<
                KEY_RSEM_EXE                  <<"=\n"<<
                KEY_GENEMARK_EXE              <<"=\n"<<
                KEY_EGGNOG_EXE                <<"=\n"<<
                KEY_EGGNOG_DOWN               <<"=\n"<<
                KEY_EGGNOG_DB                 <<"=\n"<<
                KEY_INTERPRO_EXE              <<"=\n"<<
                KEY_TAX_DB                    <<"=\n"<<
                KEY_TAX_DOWNLOAD_EXE          <<"=\n"<<
                KEY_GO_DB                     <<"=\n"<<
                KEY_GRAPH_SCRIPT              <<"=\n"
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
    if (key.compare(KEY_DIAMOND_EXE)==0)      return true;
    if (key.compare(KEY_GENEMARK_EXE)==0)     return true;
    if (key.compare(KEY_EGGNOG_EXE)==0)       return true;
    if (key.compare(KEY_EGGNOG_DOWN)==0)      return true;
    if (key.compare(KEY_EGGNOG_DB)==0)        return true;
    if (key.compare(KEY_INTERPRO_EXE)==0)     return true;
    if (key.compare(KEY_TAX_DB)==0)           return true;
    if (key.compare(KEY_GO_DB)==0)            return true;
    if (key.compare(KEY_TAX_DOWNLOAD_EXE)==0) return true;
    if (key.compare(KEY_GRAPH_SCRIPT)==0)     return true;
    return key.compare(KEY_RSEM_EXE) == 0;
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
       "EnTAP Run Information\n"   <<
       ENTAP_STATS::SOFTWARE_BREAK <<
       "Current EnTAP Version: "   << ENTAP_CONFIG::ENTAP_VERSION  <<
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
    FS_print_stats(output);
    FS_dprint(output+"\n");
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
void verify_species(boostPO::variables_map &map, SPECIES_FLAGS flag) {

    std::vector<std::string> species;
    std::string              raw_species;
    tax_serial_map_t         taxonomic_database;


    if (flag == SPECIES) {
        raw_species = map[ENTAP_CONFIG::INPUT_FLAG_SPECIES].as<std::string>();
        process_user_species(raw_species);
        species.push_back(raw_species);
    } else if (flag == CONTAMINANT) {
        species = map[ENTAP_CONFIG::INPUT_FLAG_CONTAM].as<std::vector<std::string>>();
        for (std::string &contam : species) {
            process_user_species(contam);
        }
    }
    if (species.empty()) return;

    try {
        SimilaritySearch similaritySearch = SimilaritySearch();
        taxonomic_database = similaritySearch.read_tax_map();
    } catch (const ExceptionHandler &e) {throw e;}

    for (std::string &s : species) {
        if (taxonomic_database.find(s) == taxonomic_database.end()) {
            throw ExceptionHandler("Error in one of your inputted taxons: " + s + " it is not located"
                                   " within the taxonomic database. You may remove it or select another",
                                    ENTAP_ERR::E_INPUT_PARSE);
        }
    }
    FS_dprint("Taxonomic species verified");
}



/**
 * ======================================================================
 * Description      - Finds exe paths for each piece of pipeline. WARNING: does not
 *                    check for validity of exe paths (as some parts may not want to be
 *                    ran)
 *
 * @param map       - Map of entap_config file
 * @param exe       - EnTAP execution directory
 * @return          - DIAMOND .exe path ran by enTAP::Init
 * ======================================================================
 */
void init_exe_paths(std::unordered_map<std::string, std::string> &map, std::string exe) {
    FS_dprint("Assigning execution paths. Note they are not checked for validity...");

    std::stringstream                  ss;
    std::string                        out_msg;
    boostFS::path                      exe_path(exe);
    std::pair<std::string,std::string> outpair;
    std::string temp_rsem              = map[KEY_RSEM_EXE];
    std::string temp_diamond           = map[KEY_DIAMOND_EXE];
    std::string temp_genemark          = map[KEY_GENEMARK_EXE];
    std::string temp_eggnog            = map[KEY_EGGNOG_EXE];
    std::string temp_interpro          = map[KEY_INTERPRO_EXE];
    std::string temp_eggnog_down       = map[KEY_EGGNOG_DOWN];
    std::string temp_eggnog_db         = map[KEY_EGGNOG_DB];
    std::string temp_tax_db            = map[KEY_TAX_DB];
    std::string temp_go_db             = map[KEY_GO_DB];
    std::string temp_tax_download      = map[KEY_TAX_DOWNLOAD_EXE];
    std::string temp_graphing          = map[KEY_GRAPH_SCRIPT];

    // Included software paths
    if (temp_rsem.empty())    temp_rsem           = PATHS(exe_path,Defaults::RSEM_DEFAULT_EXE);
    if (temp_diamond.empty()) temp_diamond        = PATHS(exe_path,Defaults::DIAMOND_DEFAULT_EXE);
    if (temp_genemark.empty())temp_genemark       = PATHS(exe_path,Defaults::GENEMARK_DEFAULT_EXE);
    if (temp_eggnog.empty())  temp_eggnog         = PATHS(exe_path,Defaults::EGG_EMAPPER_DEFAULT);
    if (temp_eggnog_down.empty())temp_eggnog_down = PATHS(exe_path,Defaults::EGG_DOWNLOAD_DEFAULT);
    if (temp_eggnog_db.empty())  temp_eggnog_db   = PATHS(exe_path,Defaults::EGG_SQL_DB_DEFAULT);
    if (temp_interpro.empty())   temp_interpro    = PATHS(exe_path,Defaults::INTERPRO_DEF_EXE);

    // EnTAP paths
    if (temp_tax_db.empty()) temp_tax_db     = PATHS(exe_path, ENTAP_CONFIG::TAX_DB_DEFAULT);
    if (temp_tax_download.empty()) temp_tax_download = PATHS(exe_path, Defaults::TAX_DOWNLOAD_DEF);
    if (temp_go_db.empty()) temp_go_db       = PATHS(exe_path, ENTAP_CONFIG::GO_DB_PATH_DEF);
    if (temp_graphing.empty()) temp_graphing = PATHS(exe_path, Defaults::GRAPH_SCRIPT_DEF);

    ss <<
       "\nRSEM Directory: "                  << temp_rsem         <<
       "\nGeneMarkS-T: "                     << temp_genemark     <<
       "\nDIAMOND: "                         << temp_diamond      <<
       "\nInterPro: "                        << temp_interpro     <<
       "\nEggNOG Emapper: "                  << temp_eggnog       <<
       "\nEggNOG Download: "                 << temp_eggnog_down  <<
       "\nEggNOG Database: "                 << temp_eggnog_db    <<
       "\nEnTAP Taxonomic Database: "        << temp_tax_db       <<
       "\nEnTAP Taxonomic Download Script: " << temp_tax_download <<
       "\nEnTAP Gene Ontology Database: "    << temp_go_db        <<
       "\nEnTAP Graphing Script: "           << temp_graphing;

    out_msg = ss.str();
    FS_dprint(out_msg);
    FS_print_stats(out_msg);

    DIAMOND_EXE      = temp_diamond;
    GENEMARK_EXE     = temp_genemark;
    RSEM_EXE_DIR     = temp_rsem;
    EGG_SQL_DB_PATH  = temp_eggnog_db;
    EGG_DOWNLOAD_EXE = temp_eggnog_down;
    EGG_EMAPPER_EXE  = temp_eggnog;
    INTERPRO_EXE     = temp_interpro;
    TAX_DB_PATH      = temp_tax_db;
    TAX_DOWNLOAD_EXE = temp_tax_download;
    GO_DB_PATH       = temp_go_db;
    GRAPHING_EXE     = temp_graphing;

    FS_dprint("Success! All exe paths set");
}


/**
 * ======================================================================
 * Function std::string get_exe_path(boostPO::variables_map &vm)
 *
 * Description          - Gets execution path that was used for EnTAP
 *                      - This is used for default executions with the
 *                        EnTAP config file
 *
 * Notes                - Only implemented for Unix systems now
 *
 * @param vm            - Boost map of user inputs
 * @return              - Path to executable
 * ======================================================================
 */
std::string get_exe_path(boostPO::variables_map &vm) {
    //TODO check different systems
    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_EXE_PATH)) {
        return vm[ENTAP_CONFIG::INPUT_FLAG_EXE_PATH].as<std::string>();
    }
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


/**
 * ======================================================================
 * Function bool verify_interpro(std::string database)
 *
 * Description          - Sanity check on user inputs for InterPro
 *                        databases
 *
 * Notes                - None
 *
 * @param database      - Selected databases
 * @return              - True/false is database is valid
 * ======================================================================
 */
bool verify_interpro(std::string database) {
    LOWERCASE(database);
    if (database.compare(INTER_TIGR) == 0) return true;
    if (database.compare(INTER_SFLD) == 0) return true;
    if (database.compare(INTER_PRODOM) == 0) return true;
    if (database.compare(INTER_HAMAP) == 0) return true;
    if (database.compare(INTER_PFAM) == 0) return true;
    if (database.compare(INTER_SMART) == 0) return true;
    if (database.compare(INTER_CDD) == 0) return true;
    if (database.compare(INTER_PROSITE_PROF) == 0) return true;
    if (database.compare(INTER_PROSITE_PAT) == 0) return true;
    if (database.compare(INTER_SUPERFAMILY) == 0) return true;
    if (database.compare(INTER_PRINTS) == 0) return true;
    if (database.compare(INTER_PANTHER) == 0) return true;
    if (database.compare(INTER_GENE) == 0) return true;
    if (database.compare(INTER_PIRSF) == 0) return true;
    if (database.compare(INTER_COILS) == 0) return true;
    return (database.compare(INTER_MOBI) == 0);
}


/**
 * ======================================================================
 * Function void process_user_species(std::string &input)
 *
 * Description          - Format species user has input
 *
 * Notes                - Throw error on failure
 *
 * @param input         - Species to be formatted
 * @return              - None
 * ======================================================================
 */
void process_user_species(std::string &input) {
    std::transform(input.begin(), input.end(), input.begin(), ::tolower);
    std::replace(input.begin(), input.end(), '_',' ');
}


/**
 * ======================================================================
 * Function void verify_uninformative(std::string& path)
 *
 * Description          - Sanity check on uninformative list from user
 *                      - Only checks existance/read
 *
 * Notes                - Throw error on failure
 *
 * @param path          - Path to user file
 * @return              - None
 * ======================================================================
 */
void verify_uninformative(std::string& path) {
    if (!FS_file_exists(path) || FS_file_empty(path) || !FS_file_test_open(path)) {
        throw ExceptionHandler("Path to uninformative list invalid/empty!",ENTAP_ERR::E_INPUT_PARSE);
    }
}


/**
 * ======================================================================
 * Function void verify_state(std::string &state, bool runP,
 *                            std::vector<uint16> &ontology)
 *
 * Description          - Entry to check execution paths for software based
 *                        on state
 *
 * Notes                - Throw error on failure
 *
 * @param state         - State inputted by user (or default)
 * @param runP          - Blastp flag (yes/no)
 * @param ontology      - Vector of ontology flags
 *
 * @return              - None
 * ======================================================================
 */
void verify_state(std::string &state, bool runP, std::vector<uint16> &ontology) {
    uint8 execute = 0x0;
    std::pair<bool, std::string> out;
    if (state.compare(DEFAULT_STATE) == 0) {
        execute |= DIAMOND_RUN;
        execute |= GENE_ONTOLOGY;
    }
    out = verify_software(execute, ontology);
    if (!out.first) throw ExceptionHandler(out.second, ENTAP_ERR::E_INPUT_PARSE);
}


/**
 * ======================================================================
 * Function std::pair<bool,std::string> verify_software(uint8 &states,
 *                                      std::vector<uint16> &ontology)
 *
 * Description          - Sanity check on software that will be used during
 *                        execution
 *
 * Notes                - None
 *
 * @param states        - State flags
 * @param ontology      - Vector of ontology flags
 *
 * @return              - Pair of yes/no failure and error msg string
 * ======================================================================
 */
std::pair<bool,std::string> verify_software(uint8 &states,std::vector<uint16> &ontology) {
    FS_dprint("Verifying software...");

    if (states & DIAMOND_RUN) {
        if (!FS_file_exists(TAX_DB_PATH) || FS_file_empty(TAX_DB_PATH))
            return std::make_pair(false, "Could not find the taxonomic database");
        if (!SimilaritySearch::is_executable()) {
            return std::make_pair(false, "Could not execute a test run of DIAMOND, be sure"
                    " it's properly installed and the path is correct");
        }
    }
    if (states & GENE_ONTOLOGY) {
        if (!FS_file_exists(GO_DB_PATH) || FS_file_empty(GO_DB_PATH))
            return std::make_pair(false, "Could not find Gene Ontology database or invalid");
        for (uint16 flag : ontology) {
            switch (flag) {
                case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
                    if (!FS_file_exists(EGG_SQL_DB_PATH))
                        return std::make_pair(false, "Could not find EggNOG SQL database");
                    if (!FS_file_exists(EGG_EMAPPER_EXE) || !ModEggnog::is_executable())
                        return std::make_pair(false, "Could not find or test EggNOG Emapper, "
                                "ensure python is properly installed and the paths are correct");
                    break;
                case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
                    // TODO
                    break;
                default:
                    break;
            }
        }
    }

    FS_dprint("Success!");
    return std::make_pair(true, "");
}