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
#include "common.h"
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
#include <thread>
#include "QuerySequence.h"
#include "FrameSelection.h"
#include "ExpressionAnalysis.h"
#include "SimilaritySearch.h"
#include "FileSystem.h"
#include <boost/regex.hpp>
#include <queue>
#include <iomanip>
#include "UserInput.h"
//**************************************************************


namespace entapExecute {

    //*********************** Local Variables ***********************
    ExecuteStates           executeStates;
    std::string             _outpath;
    std::string             _entap_outpath;
    bool                    _blastp;          // false for blastx, true for _blastp
    std::string             _input_path;      // FASTA changes depending on execution
    int                     _threads;
    databases_t             _databases;       // NCBI+UNIPROT+Other
    FileSystem             *_pFileSystem;
    UserInput              *_pUserInput;
    //**************************************************************

    //******************** Local Prototype Functions ***************
    std::string filter_transcriptome(std::string &);
    void verify_state(std::queue<char> &, bool &);
    bool valid_state(enum ExecuteStates);
    void exit_error(ExecuteStates);
    std::vector<std::string> verify_databases(std::vector<std::string>, std::vector<std::string>,
                                              std::vector<std::string>, std::string);

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
    void execute_main(UserInput *user_input, FileSystem *filesystem) {
        FS_dprint("EnTAP Executing...");

        std::vector<uint16>                     ontology_flags;
        std::string                             original_input;  // ALWAYS use for Expression
        std::queue<char>                        state_queue;
        bool                                    state_flag;

        if (user_input == nullptr || filesystem == nullptr) {
            throw ExceptionHandler("Unable to allocate memory to EnTAP Execution", ERR_ENTAP_INPUT_PARSE);
        }

        executeStates           = INIT;
        state_flag              = false;
        _pUserInput             = user_input;
        _pFileSystem            = filesystem;

        // Pull relevant info input by the user
        _input_path    = _pUserInput->get_user_input<std::string>(UInput::INPUT_FLAG_TRANSCRIPTOME);
        original_input = _input_path;
        _blastp        = _pUserInput->has_input(UInput::INPUT_FLAG_RUNPROTEIN);
        ontology_flags = _pUserInput->get_user_input<std::vector<uint16>>(UInput::INPUT_FLAG_ONTOLOGY);
        state_queue    = _pUserInput->get_state_queue();    // Will NOT be empty
        _databases     = _pUserInput->get_user_input<databases_t>(UInput::INPUT_FLAG_DATABASE);

        // Set/create outpaths
        _outpath       = _pFileSystem->get_root_path();
        _entap_outpath = PATHS(_outpath, ENTAP_OUTPUT);     // transcriptome outpath, original will be copied
        _pFileSystem->create_dir(_entap_outpath);
        _pFileSystem->create_dir(_outpath);

        try {
            verify_state(state_queue, state_flag);         // Set state transition
            // Read input transcriptome
            QueryData *pQUERY_DATA = new QueryData(
                    _input_path,        // User transcriptome
                    _entap_outpath,     // Transcriptome directory
                    _pUserInput,        // User input map
                    _pFileSystem);      // Filesystem object
            GraphingManager graphingManager = GraphingManager(GRAPHING_EXE);

            while (executeStates != EXIT) {
                switch (executeStates) {
                    case FRAME_SELECTION: {
                        FS_dprint("STATE - FRAME SELECTION");
                        std::unique_ptr<FrameSelection> frame_selection(new FrameSelection(
                                _input_path,
                                _pFileSystem,
                                _pUserInput,
                                &graphingManager,
                                pQUERY_DATA
                        ));
                        if ((_blastp && pQUERY_DATA->DATA_FLAG_GET(QueryData::IS_PROTEIN))) {
                            FS_dprint("Protein sequences input, skipping frame selection");
                        } else if (!_blastp) {
                            FS_dprint("Blastx selected, skipping frame selection");
                        } else {
                            _input_path = frame_selection->execute(_input_path);
                            pQUERY_DATA->DATA_FLAG_SET(QueryData::IS_PROTEIN);
                            pQUERY_DATA->DATA_FLAG_SET(QueryData::SUCCESS_FRAME_SEL);
                        }
                    }
                        break;
                    case EXPRESSION_FILTERING: {
                        FS_dprint("STATE - EXPRESSION FILTERING");
                        if (!_pUserInput->has_input(UInput::INPUT_FLAG_ALIGN)) {
                            FS_dprint("No alignment file specified, skipping expression analysis");
                        } else {
                            std::unique_ptr<ExpressionAnalysis> expression(new ExpressionAnalysis(
                                original_input,
                                &graphingManager,
                                pQUERY_DATA,
                                _pFileSystem,
                                _pUserInput
                            ));
                            _input_path = expression->execute(original_input);
                            pQUERY_DATA->DATA_FLAG_SET(QueryData::SUCCESS_EXPRESSION);
                        }
                    }
                        break;
                    case FILTER:
                        _input_path = filter_transcriptome(_input_path);  // Just copies final transcriptome
                        break;
                    case SIMILARITY_SEARCH: {
                        FS_dprint("STATE - SIM SEARCH RUN");
                        // Spawn sim search object
                        std::unique_ptr<SimilaritySearch> sim_search(new SimilaritySearch(
                                _databases,
                                _input_path,
                                _pUserInput,
                                _pFileSystem,
                                &graphingManager,
                                pQUERY_DATA
                        ));

                        sim_search->execute(_input_path, _blastp);
                        pQUERY_DATA->DATA_FLAG_SET(QueryData::SUCCESS_SIM_SEARCH);
                        FS_dprint("STATE - SIM SEARCH PARSE");
                        sim_search->parse_files(_input_path);
                        break;
                    }
                    case GENE_ONTOLOGY: {
                        FS_dprint("STATE - GENE ONTOLOGY");
                        std::unique_ptr<Ontology> ontology(new Ontology(
                                _input_path,
                                _pUserInput,
                                &graphingManager,
                                pQUERY_DATA,
                                _pFileSystem
                        ));
                        ontology->execute(_input_path);
                        pQUERY_DATA->DATA_FLAG_SET(QueryData::SUCCESS_ONTOLOGY);
                        break;
                    }
                    default:
                        executeStates = EXIT;
                        break;
                }
                verify_state(state_queue, state_flag);
            }

            // *************************** Exit Stuff ********************** //
            pQUERY_DATA->final_statistics(_outpath, ontology_flags);
            _pFileSystem->directory_iterate(true, _outpath);   // Delete empty files
            SAFE_DELETE(pQUERY_DATA);
        } catch (const ExceptionHandler &e) {
            exit_error(executeStates);
            throw e;
        }
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
        FS_dprint("Beginning to copy final transcriptome to be used...");

        boostFS::path file_name;
        std::string   file_name_str;
        std::string   out_path;

        file_name = input_path;
        file_name_str = file_name.filename().stem().string() + FINAL_TRANSCRIPTOME_TAG;
        out_path = PATHS(_entap_outpath, file_name_str);
        boostFS::copy_file(input_path,out_path,boostFS::copy_option::overwrite_if_exists);

        FS_dprint("Success!");
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
        FS_dprint("verifying state...");
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
        FS_dprint("Success!");
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
        return (s >= EXPRESSION_FILTERING && s <= EXIT);
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
#if 0
        switch (s) {
            case INIT:
                break;
            case EXPRESSION_FILTERING:
                ss <<
                   "EnTAP failed execution during the Expression Filtering stage with\n"
                           "previous stages executing correctly.\n";
                break;
            case FRAME_SELECTION:
                ss <<
                   "EnTAP failed execution during the  Frame Selection stage with\n"
                           "previous stages executing correctly.\n";
                break;
            case FILTER:
                break;
            case SIMILARITY_SEARCH:
                ss <<
                   "EnTAP failed execution during the Similarity Searching \n"
                           " with previous stages executing correctly.\n"
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
#endif
        ss <<
           "Here are a few ways to help diagnose some general issues:\n"
                   "\t1. Check the (detailed) printed error message below\n"
                   "\t2. Review the .err files of the execution stage (they will\n"
                   "\t\tbe in the directory for whatever stage you failed on\n"
                   "\t3. Check the debug.txt file that is printed after execution\n"
                   "\t4. Ensure your paths/inputs are correct in entap_config.txt\n"
                   "\t\tand log_file.txt (this will show your inputs)\n";
        ss << "------------------------------------";
        std::cerr<<ss.str()<<std::endl;
    }
}
