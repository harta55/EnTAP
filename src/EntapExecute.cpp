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
#include <boost/serialization/unordered_map.hpp>
#include <fstream>
#include <map>
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
#include <thread>
#include "QuerySequence.h"
#include "FrameSelection.h"
#include "ExpressionAnalysis.h"
#include "SimilaritySearch.h"
#include <boost/regex.hpp>
#include <queue>
#include <iomanip>
//**************************************************************


namespace entapExecute {

    //*********************** Globals *****************************

    ExecuteStates           executeStates;
    std::string             _outpath;
    std::string             _entap_outpath;
    bool                    _EXPRESSION_SUCCESS;     // True if this stage was ran/success
    bool                    _FRAME_SELETION_SUCCESS;
    bool                    _SIM_SEARCH_SUCCESS;
    bool                    _ONTOLOGY_SUCCESS;
    bool                    _blastp;          // false for blastx, true for _blastp
    std::string             _input_path;      // FASTA changes depending on execution
    int                     _threads;
    std::vector<std::string>_databases;       // NCBI+UNIPROT+Other

    //**************************************************************

/**
 * ======================================================================
 * Function execute_main(boostPO::variables_map      &user_input,
 *                       std::string                 exe_path,
                         std::unordered_map<std::string, std::string> &config_map)
 *
 * Description          - Entry into main annotation portion of EnTAP, manages
 *                        calling each state of the pipeline
 *                      - Calculates overall statistics and parses input
 *                        transcriptome
 *                      - Parses input databases
 *                      - Determines execution paths based on configuration file
 *
 * Notes                - Entry
 *
 * @param user_input    - Boost parsed user input flags
 * @param exe_path      - Path to EnTAP executable and main directory
 * @param config_map    - Map of the EnTAP configuration file
 *
 * @return              - None
 *
 * =====================================================================
 */
    void execute_main(boost::program_options::variables_map &user_input) {
        print_debug("EnTAP Executing...");

        std::vector<uint16>                     ontology_flags;
        std::vector<std::string>                other_databases; // -d Command databases
        std::pair<std::string,std::string>      diamond_pair;    // best_hits.fa,no_hits.fa
        std::string                             no_database_hits;// No DIAMOND
        std::string                             original_input;  // ALWAYS use for Expression
        std::queue<char>                        state_queue;
        bool                                    trim_flag;       // User trim flag selected
        bool                                    state_flag;
        bool                                    is_complete;     // All input sequences are complete genes

        executeStates           = INIT;
        state_flag              = false;
        _EXPRESSION_SUCCESS     = false;
        _FRAME_SELETION_SUCCESS = false;
        _SIM_SEARCH_SUCCESS     = false;
        _ONTOLOGY_SUCCESS       = false;

        _input_path    = user_input[ENTAP_CONFIG::INPUT_FLAG_TRANSCRIPTOME].as<std::string>();
        _threads       = get_supported_threads(user_input);
        _blastp        = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_RUNPROTEIN);
        original_input = _input_path;
        trim_flag      = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_TRIM);
        is_complete    = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_COMPLETE);
        diamond_pair   = std::make_pair(_input_path,"");
        ontology_flags = user_input[ENTAP_CONFIG::INPUT_FLAG_ONTOLOGY].as<std::vector<uint16>>();

        boostFS::path working_dir(boostFS::current_path());
        _outpath       = PATHS(working_dir, user_input["tag"].as<std::string>());
        _entap_outpath = PATHS(_outpath, ENTAP_OUTPUT);
        boostFS::create_directories(_entap_outpath);
        boostFS::create_directories(_outpath);

        // init databases
        if (user_input.count("database")) {
            other_databases = user_input["database"].as<std::vector<std::string>>();
        } else other_databases.push_back(ENTAP_CONFIG::NCBI_NULL);

        // init state control
        if (user_input.count("state")) {
            std::string user_state_str = user_input["state"].as<std::string>();
            for (char c : user_state_str) {
                state_queue.push(c);
            }
        }

        try {
            _databases = verify_databases(user_input["uniprot"].as<std::vector<std::string>>(),
                                         user_input["ncbi"].as<std::vector<std::string>>(),
                                         other_databases, "");
            verify_state(state_queue, state_flag);
            QueryData QUERY_DATA = QueryData(_input_path, _entap_outpath,is_complete, trim_flag);
            GraphingManager graphingManager = GraphingManager(GRAPHING_EXE);
            std::unique_ptr<SimilaritySearch> sim_search(new SimilaritySearch(
                    _databases, _input_path, _threads, _outpath,
                    user_input, &graphingManager, &QUERY_DATA
            ));

            while (executeStates != EXIT) {
                switch (executeStates) {
                    case FRAME_SELECTION: {
                        print_debug("STATE - FRAME SELECTION");
                        std::unique_ptr<FrameSelection> frame_selection(new FrameSelection(
                                _input_path, _outpath, user_input, &graphingManager, &QUERY_DATA
                        ));
                        if ((_blastp && QUERY_DATA.is_protein())) {
                            print_debug("Protein sequences input, skipping frame selection");
                            QUERY_DATA.flag_transcripts(FRAME_SELECTION);
                        } else if (!_blastp) {
                            print_debug("Blastx selected, skipping frame selection");
                            QUERY_DATA.flag_transcripts(FRAME_SELECTION);
                        } else {
                            _input_path = frame_selection->execute(_input_path);
                            _FRAME_SELETION_SUCCESS = true;
                        }
                        frame_selection.release();
                    }
                        break;
                    case RSEM: {
                        print_debug("STATE - EXPRESSION");
                        std::unique_ptr<ExpressionAnalysis> expression(new ExpressionAnalysis(
                                original_input, _threads, _outpath, user_input, &graphingManager, &QUERY_DATA
                        ));
                        if (!user_input.count(ENTAP_CONFIG::INPUT_FLAG_ALIGN)) {
                            print_debug("No alignment file specified, skipping expression analysis");
                            QUERY_DATA.flag_transcripts(RSEM);
                        } else {
                            _input_path = expression->execute(original_input);
                            _EXPRESSION_SUCCESS = true;
                        }
                        expression.release();
                    }
                        break;
                    case FILTER:
                        _input_path = filter_transcriptome(_input_path);
                        break;
                    case DIAMOND_RUN:
                        print_debug("STATE - SIM SEARCH RUN");
                        sim_search->execute(_input_path, _blastp);
                        _SIM_SEARCH_SUCCESS = true;
                        break;
                    case DIAMOND_PARSE:
                        print_debug("STATE - SIM SEARCH PARSE");
                        diamond_pair = sim_search->parse_files(_input_path);
                        _input_path = diamond_pair.first;
                        no_database_hits = diamond_pair.second;
                        break;
                    case GENE_ONTOLOGY: {
                        print_debug("STATE - GENE ONTOLOGY");
                        std::unique_ptr<Ontology> ontology(new Ontology(
                                _threads, _outpath, _input_path, user_input, &graphingManager,
                                &QUERY_DATA, _blastp
                        ));
                        ontology->execute(_input_path, no_database_hits);
                        _ONTOLOGY_SUCCESS = true;
                        ontology.release();
                    }
                        break;
                    default:
                        executeStates = EXIT;
                        break;
                }
                verify_state(state_queue, state_flag);
            }
            QUERY_DATA.set_EXPRESSION_SUCCESS(_EXPRESSION_SUCCESS);
            QUERY_DATA.set_FRAME_SELECTION_SUCCESS(_FRAME_SELETION_SUCCESS);
            QUERY_DATA.set_ONTOLOGY_SUCCESS(_ONTOLOGY_SUCCESS);
            QUERY_DATA.set_SIM_SEARCH_SUCCESS(_SIM_SEARCH_SUCCESS);
            QUERY_DATA.final_statistics(_outpath, ontology_flags);
        } catch (const ExceptionHandler &e) {
            exit_error(executeStates);
            throw e;
        }
    }


    std::vector<std::string> verify_databases(std::vector<std::string> uniprot, std::vector<std::string> ncbi,
                                            std::vector<std::string> database, std::string exe) {
        print_debug("Verifying databases...");
        // return file paths
        // config file paths already exist (checked in main)
        std::vector<std::string>        file_paths;
        std::string                     path;
        std::string                     config_path;

#if NCBI_UNIPROT
        print_debug("Verifying uniprot databases...");
        if (uniprot.size() > 0) {
            for (auto const &u_flag:uniprot) {
                if (u_flag.compare(ENTAP_CONFIG::INPUT_UNIPROT_NULL) != 0) {
                    if (u_flag == ENTAP_CONFIG::INPUT_UNIPROT_SWISS) {
                        config_path = config.at(ENTAP_CONFIG::KEY_UNIPROT_SWISS);
                    } else if (u_flag == ENTAP_CONFIG::INPUT_UNIPROT_TREMBL) {
                        config_path = config.at(ENTAP_CONFIG::KEY_UNIPROT_TREMBL);
                    } else if (u_flag == ENTAP_CONFIG::INPUT_UNIPROT_UR90) {
                        config_path = config.at(ENTAP_CONFIG::KEY_UNIPROT_UR90);
                    } else if (u_flag == ENTAP_CONFIG::INPUT_UNIPROT_UR100) {
                        config_path = config.at(ENTAP_CONFIG::KEY_UNIPROT_UR100);
                    }
                    if (!config_path.empty()) {
                        print_debug("Config file database found, using this path at: " +
                                             config_path);
                        path = config_path;
                    } else {
                        path = exe + ENTAP_CONFIG::UNIPROT_INDEX_PATH + u_flag + ".dmnd";
                    }
                    if (!file_exists(path))
                        throw ExceptionHandler("Database located at: " + path + " not found", ENTAP_ERR::E_INPUT_PARSE);
                    file_paths.push_back(path);
                } else {
                    print_debug("No/null Uniprot databases detected");
                    break;
                }
            }
        }
        print_debug("Complete");
        print_debug("Verifying NCBI databases...");
        if (ncbi.size() > 0) {
            for (auto const &u_flag:ncbi) {
                if (u_flag.compare(ENTAP_CONFIG::NCBI_NULL) != 0) {
                    if (u_flag == ENTAP_CONFIG::NCBI_NONREDUNDANT) {
                        config_path = config.at(ENTAP_CONFIG::KEY_NCBI_NR);
                    } else if (u_flag == ENTAP_CONFIG::NCBI_REFSEQ_PLANT) {
                        config_path = config.at(ENTAP_CONFIG::KEY_NCBI_REFSEQ_SEPARATE);
                    } else if (u_flag == ENTAP_CONFIG::NCBI_REFSEQ_COMP) {
                        config_path = config.at(ENTAP_CONFIG::KEY_NCBI_REFSEQ_COMPLETE);
                    }
                    if (!config_path.empty()) {
                        print_debug("Config file database found, using this path at: " +
                                             config_path);
                        path = config_path;
                    } else {
                        path = exe + ENTAP_CONFIG::NCBI_INDEX_PATH + u_flag + ".dmnd";
                    }
                    if (!file_exists(path))
                        throw ExceptionHandler("Database located at: " + path + " not found", ENTAP_ERR::E_INPUT_PARSE);
                    file_paths.push_back(path);
                } else {
                    print_debug("No/null NCBI databases detected");
                    break;
                }
            }
        }
        print_debug("Complete");

#endif
        print_debug("Verifying other databases...");
        if (database.size() > 0) {
            for (auto const &data_path:database) {
                if (data_path.compare(ENTAP_CONFIG::NCBI_NULL) == 0) continue;
                if (!file_exists(data_path)) {
                    throw ExceptionHandler("Database located at: " + data_path + " not found",
                                           ENTAP_ERR::E_INPUT_PARSE);
                }
                boostFS::path bpath(data_path);
                std::string ext = bpath.extension().string();
                if (ext.compare(".dmnd") == 0) {
                    print_debug("User has input a diamond indexed database at: " + data_path);
                    file_paths.push_back(data_path);
                    continue;
                } else {
                    //todo fix not really used yet
                    print_debug("User has input a database at: " + data_path);
                    std::string test_path = PATHS(exe,ENTAP_CONFIG::BIN_PATH) + data_path + ".dmnd";
                    print_debug("Checking if indexed file exists at: " + test_path);
                    if (!file_exists(test_path)) {
                        throw ExceptionHandler("Database located at: " + data_path + " not found",
                                               ENTAP_ERR::E_INPUT_PARSE);
                    } else {
                        file_paths.push_back(test_path);
                    }
                }
            }
        }
        print_debug("Verification complete!");
        if (file_paths.size() > 0) {
            std::string database_final = "\n\nDatabases selected:\n";
            for (std::string base: file_paths) {
                database_final += base + "\n";
            }
            print_debug(database_final);
        } else {
            print_debug("No databases selected, some functionality may not be able to run");
        }
        return file_paths;
    }


    /**
     * ======================================================================
     * Function std::string filter_transcriptome(std::string    &input_path)
     *
     * Description          - Merely selects transcriptome that will
     *                        continue in pipeline and copies it to entap_out directory
     *
     * Notes                - None
     *
     * @param input_path    - Input transcriptome (expression and/or frame selected)
     * @return              - Copied transcriptome
     *
     * =====================================================================
     */
    std::string filter_transcriptome(std::string &input_path) {
        print_debug("Beginning to copy final transcriptome to be used...");

        boostFS::path file_name;
        std::string   file_name_str;
        std::string   out_path;

        file_name = input_path;
        file_name_str = file_name.filename().stem().string() + "_final.fasta";
        out_path = (boostFS::path(_entap_outpath) / file_name_str).string();
        boostFS::copy_file(input_path,out_path,boostFS::copy_option::overwrite_if_exists);

        print_debug("Success!");
        return out_path;
    }




/**
 * ======================================================================
 * Function verify_state(std::queue<char> &queue, bool &test)
 *
 * Description          - Computes the next state that will be executed
 *                        based upon --state flag (default: +)
 *                      - User input not thoroughly checked and is experimental
 *                      - Normal functionality just transitions between
 *                        states as normal
 *
 * Notes                - Only assuming between 0-9, no digit states
 *
 * @param queue         - State queue as inputted by user (or default '+')
 * @param test          - Flag used if previous character was a '+'
 * @return              - None
 * ======================================================================
 */
    void verify_state(std::queue<char> &queue, bool &test) {
        print_debug("verifying state...");
        if (queue.empty()) {
            executeStates = static_cast<ExecuteStates>(executeStates + 1);
            if (!valid_state(executeStates)) executeStates = EXIT;
            return;
        }
        char first = queue.front();
        if (first == 'x') {
            queue.pop();
            test = false;
            if (queue.empty()) {
                executeStates = EXIT;
                return;
            }
            verify_state(queue, test);
            return;
        }
        if (first == '+') {
            // assuming proper cast, state has been evaluated before
            test = true;
            queue.pop();
            char second = queue.front(); // assuming number
            if (!second) {
                // end of queue, might handle differently
                executeStates = static_cast<ExecuteStates>(executeStates + 1);
                if (!valid_state(executeStates)) executeStates = EXIT;
                return;
            }
            verify_state(queue, test);
        } else {
            // some number (assuming it was parsed before)
            int i = first - '0';
            if (executeStates < i && test) {
                executeStates = static_cast<ExecuteStates>(executeStates + 1);
                return;
            } else if (executeStates < i && !test) {
                executeStates = static_cast<ExecuteStates>(i);
                return;
            }
            if (executeStates == i) {
                queue.pop();
                test = false;
                verify_state(queue, test);
            }
        }
        print_debug("Success!");
    }


/**
 * ======================================================================
 * Function bool valid_state(ExecuteStates s)
 *
 * Description          - Ensures computed state from verify_state() func
 *                        is valid and within the range
 *
 * Notes                - None
 *
 * @param s             - State to verify
 * @return              - None
 * ======================================================================
 */
    bool valid_state(ExecuteStates s) {
        return (s >= RSEM && s <= EXIT);
    }


    /**
     * ======================================================================
     * Function void exit_error(ExecuteStates s)
     *
     * Description          - Prints a "summary" report for user during a
     *                        fatal error depending on stage
     *
     * Notes                - None
     *
     * @param s             - Execution state
     * @return              - None
     * ======================================================================
     */
    void exit_error(ExecuteStates s) {
        std::stringstream ss;

        ss << "------------------------------------\n";
        switch (s) {
            case INIT:
                break;
            case RSEM:
                ss <<
                   "EnTAP failed execution during the Expression Filtering stage with\n"
                           "previous stages executing correctly.\n";
                break;
            case FRAME_SELECTION:
                ss <<
                   "EnTAP failed execution during the Expression Filtering stage with\n"
                           "previous stages executing correctly.\n";
                break;
            case FILTER:
                break;
            case DIAMOND_RUN:
                ss <<
                   "EnTAP failed execution during executing Similarity Searching with\n"
                           "previous stages executing correctly.\n";
                break;
            case DIAMOND_PARSE:
                ss <<
                   "EnTAP failed execution during the parsing of Similarity Searching \n"
                           "data with previous stages executing correctly.\n"
                           "It seems that the similarity searching worked, but \n"
                           "some issue in parsing caused an error!\n";
                break;
            case GENE_ONTOLOGY:
                ss <<
                   "EnTAP failed execution during the Gene Ontology stage with\n"
                           "previous stages executing correctly.\n";
                break;
            default:
                break;
        }
        ss <<
           "Here are a few ways to help diagnose some general issues:\n"
                   "\t1. Check the (detailed) printed error message below\n"
                   "\t2. Review the .err files of the execution stage (they will\n"
                   "\t\tbe in the directory for whateve stage you failed on\n"
                   "\t3. Check the debug.txt file that is printed after execution\n"
                   "\t4. Ensure your paths/inputs are correct in entap_config.txt\n"
                   "\t\tand log_file.txt (this will show your inputs)\n";
        ss << "------------------------------------";
        std::cerr<<ss.str()<<std::endl;
    }
}
