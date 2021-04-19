/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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
#include "EntapExecute.h"
//**************************************************************


namespace entapExecute {

    //*********************** Local Variables ***********************
    ExecuteStates           executeStates;    // Current Execution state
    FileSystem             *pFileSystem;      // Pointer to EnTAP filesystem
    UserInput              *pUserInput;       // Pointer to input flags from user
    //**************************************************************

    //******************** Local Prototype Functions ***************
    void verify_state(std::queue<char> &queue, bool test);
    bool valid_state(enum ExecuteStates current_state);
    void exit_error(ExecuteStates exiting_state);
    //**************************************************************

/**
 * ======================================================================
 * Function execute_main(UserInput *user_input, FileSystem *filesystem)
 *
 * Description          - Entry into main annotation portion of EnTAP, manages
 *                        calling each state of the pipeline
 *                      - Calculates overall statistics and parses input
 *                        transcriptome
 *
 * Notes                - Entry
 *
 * @param user_input    - Pointer to user input flags
 * @param filesystem    - Pointer to EnTAP filesystem
 *
 * @return              - None
 *
 * =====================================================================
 */
    void execute_main(UserInput *user_input, FileSystem *filesystem) {
        FS_dprint("EnTAP Executing...");

        std::string                             original_input;  // Absolute path to ORIGINAL transcriptome (ALWAYS use for Expression)
        std::string                             input_path;      // Absolute path to transciptome (WARNING changes depending on user selection)
        std::string                             final_out_dir;   // Absolute path to output directory for final stats
        std::string                             transcriptome_outpath;
        std::string                             transcriptome_dir;
        std::string                             transcrtipeom_out_filename;
        std::string                             entap_graphing_path;
        ent_input_multi_str_t                   databases;       // NCBI+UNIPROT+Other
        std::queue<char>                        state_queue;
        bool                                    state_flag;
        bool                                    blastp;                 // false for blastx, true for mBlastp
        vect_uint16_t                           entap_database_types;   // Database types from user
        EntapDatabase::DATABASE_TYPE            entap_database_type;    // First database type selected by user that will be used
        std::string                             entap_database_path;    // Path to EnTAP database we are using
        EntapDataPtrs                           entap_data_ptrs;        // Struct of Execution pointers
        EntapDatabase*                          pEntap_Database=nullptr; // Pointer to the EnTAP database (binary / sql)
        QueryData*                              pQUERY_DATA=nullptr;    // Pointer to the trancsriptome data
        GraphingManager*                        pGraphing_Manager=nullptr;   // Pointer to the graphing manager

        // Ensure data exists
        if (user_input == nullptr || filesystem == nullptr) {
            throw ExceptionHandler("Unable to allocate memory to EnTAP Execution", ERR_ENTAP_MEM_ALLOC);
        }

        executeStates           = INIT;
        state_flag              = false;
        pUserInput              = user_input;
        pFileSystem             = filesystem;
        entap_data_ptrs         = EntapDataPtrs();

        // Pull relevant info input by the user
        input_path     = pUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_TRANSCRIPTOME);
        original_input = input_path;
        blastp         = pUserInput->has_input(INPUT_FLAG_RUNPROTEIN);
        state_queue    = pUserInput->get_state_queue();    // Will NOT be empty, default is +
        databases     = pUserInput->get_user_input<ent_input_multi_str_t>(INPUT_FLAG_DATABASE);
        entap_graphing_path = pUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_ENTAP_GRAPH);

        // Find database type that will be used by the rest (use 0 index no matter what)
        entap_database_types = pUserInput->get_user_input<ent_input_multi_int_t>(INPUT_FLAG_DATABASE_TYPE);
        entap_database_type  = static_cast<EntapDatabase::DATABASE_TYPE>(entap_database_types[0]);
        entap_database_path  = pUserInput->get_entap_database_path(entap_database_type);

        // Set/create outpaths
        final_out_dir  = pFileSystem->get_final_outdir();
        pFileSystem->create_transcriptome_dir();        // Directory where filtered transcriptomes will be copied
        transcrtipeom_out_filename = pFileSystem->get_filename(input_path, true);   // Pull filename from input transcriptome
        transcriptome_dir     = pFileSystem->get_trancriptome_dir();
        transcriptome_outpath = PATHS(transcriptome_dir, transcrtipeom_out_filename);

        try {
            verify_state(state_queue, state_flag);         // Set state transition

            // Initialize Query Data
            pQUERY_DATA = new QueryData(
                    input_path,        // User transcriptome
                    transcriptome_outpath,     // Transcriptome directory
                    pUserInput,        // User input map
                    pFileSystem);      // Filesystem object

            // Initialize Graphing Manager
            pGraphing_Manager = new GraphingManager(entap_graphing_path, pFileSystem);

            // Initialize EnTAP database
            pEntap_Database = new EntapDatabase(filesystem);
            if (!pEntap_Database->set_database(entap_database_type, entap_database_path)) {
                throw ExceptionHandler("Unable to initialize EnTAP database\n" +
                                       pEntap_Database->print_error_log(), ERR_ENTAP_READ_ENTAP_DATA_GENERIC);
            }

            // Compile all data pointers needed during execution
            entap_data_ptrs.mpEntapDatabase = pEntap_Database;
            entap_data_ptrs.mpFileSystem   = filesystem;
            entap_data_ptrs.mpUserInput    = user_input;
            entap_data_ptrs.mpGraphingManager = pGraphing_Manager;
            entap_data_ptrs.mpQueryData    = pQUERY_DATA;

            if (entap_data_ptrs.is_null()) {
                throw ExceptionHandler("Unable to allocate memory", ERR_ENTAP_MEM_ALLOC);
            }

            while (executeStates != EXIT) {
                switch (executeStates) {

                    case EXPRESSION_FILTERING: {
                        FS_dprint("STATE - EXPRESSION FILTERING");
                        if (!user_input->run_expression_filtering()) {
                            FS_dprint("No alignment file specified, skipping expression analysis");
                            pQUERY_DATA->header_set(ENTAP_HEADER_EXP_FPKM, false);
                            pQUERY_DATA->header_set(ENTAP_HEADER_EXP_TPM, false);
                        } else {
                            FS_dprint("Continuing with Expression Analysis");
                            // Proceed with expression analysis
                            std::unique_ptr<ExpressionAnalysis> expression(new ExpressionAnalysis(
                                    original_input, entap_data_ptrs
                            ));
                            // WARNING: set out next input path to the filtered one
                            input_path = expression->execute(original_input);

                            // Copy filtered file to entap transcriptome directory
                            std::string transc_filter_filename = pUserInput->get_user_transc_basename() + TRANSCRIPTOME_FILTERED_TAG;
                            std::string transc_filter_outpath  = PATHS(transcriptome_dir, transc_filter_filename);

                            pFileSystem->delete_file(transc_filter_outpath);
                            uint32 sequence_flags = 0;
                            QueryData::SEQUENCE_TYPES sequence_type;
                            sequence_flags |= QuerySequence::QUERY_EXPRESSION_KEPT;
                            pQUERY_DATA->is_protein_data() ? sequence_type = QueryData::SEQUENCE_AMINO_ACID : sequence_type = QueryData::SEQUENCE_NUCLEOTIDE;
                            if (pQUERY_DATA->print_transcriptome(sequence_flags, transc_filter_outpath, sequence_type)){
                                FS_dprint("Expression filtered transcriptome generated to: " + transc_filter_outpath);
                                // WARNING: Set our next input path to the frame selected version
                                input_path = transc_filter_outpath;
                            } else {
                                throw ExceptionHandler("ERROR: unable to generate filtered transcriptome from expression analysis results",
                                                       ERR_ENTAP_GENERATE_TRANSCRIPTOME);
                            }
                        }
                        break;
                    }

                    case FRAME_SELECTION: {
                        FS_dprint("STATE - FRAME SELECTION");
                        bool run_frame_selection;
                        if (pUserInput->run_frame_selection(pQUERY_DATA, run_frame_selection)) {
                            // We were able to determine if we want to run frame selection
                            if (run_frame_selection) {

                                FS_dprint("Continuing with frame selection process...");
                                std::unique_ptr<FrameSelection> frame_selection(new FrameSelection(
                                        input_path, entap_data_ptrs
                                ));

                                frame_selection->execute(input_path);

                                // Copy frame selected file to the trancriptome directory
                                std::string transc_protein_filename = pUserInput->get_user_transc_basename() + TRANSCRIPTOME_FRAME_TAG;
                                std::string transc_protein_outpath  = PATHS(transcriptome_dir, transc_protein_filename);
                                uint32 sequence_flags=0;
                                sequence_flags |= QuerySequence::QUERY_FRAME_KEPT;
                                pFileSystem->delete_file(transc_protein_outpath);
                                if (pQUERY_DATA->print_transcriptome(sequence_flags, transc_protein_outpath,
                                                                     QueryData::SEQUENCE_AMINO_ACID)){
                                    FS_dprint("Protein transcriptome generated to: " + transc_protein_outpath);
                                    // WARNING: Set our next input path to the frame selected version
                                    input_path = transc_protein_outpath;
                                } else {
                                    throw ExceptionHandler("ERROR: unable to generate protein transcriptome from frame selecttion results",
                                                           ERR_ENTAP_GENERATE_TRANSCRIPTOME);
                                }
                            } else {
                                // No, we are going to skip frame selection process
                                pQUERY_DATA->header_set(ENTAP_HEADER_FRAME, false);
                            }
                        } else {
                            // Should never happen
                            throw ExceptionHandler("Unable to determine if we want to run Frame Selection",
                                ERR_ENTAP_INPUT_PARSE);
                        }

                        break;
                    }

                    case COPY_FINAL_TRANSCRIPTOME: {
                        // At this point, input_path is the final transcriptome we will be using
                        // copy this final trancriptome and then set input_path to the copied trancriptome
                        // absolute paths only!
                        FS_dprint("Beginning to copy final transcriptome to be used...");

                        std::string file_name = pUserInput->get_user_transc_basename() + TRANSCRIPTOME_FINAL_TAG;
                        std::string out_path = PATHS(transcriptome_dir, file_name);
                        uint32 sequence_flags = 0;
                        QueryData::SEQUENCE_TYPES sequence_type;
                        sequence_flags |= QuerySequence::QUERY_FRAME_KEPT;
                        sequence_flags |= QuerySequence::QUERY_EXPRESSION_KEPT;
                        blastp ? sequence_type = QueryData::SEQUENCE_AMINO_ACID
                               : sequence_type = QueryData::SEQUENCE_NUCLEOTIDE;
                       pFileSystem->delete_file(out_path);
                        if (pQUERY_DATA->print_transcriptome(sequence_flags, out_path, sequence_type)) {
                            FS_dprint("FINAL transcriptome generated to: " + out_path);
                            // WARNING: Set our next input path to the frame selected version
                            input_path = out_path;
                        } else {
                            throw ExceptionHandler(
                                    "ERROR: unable to generate final transcriptome from frame selecttion results",
                                    ERR_ENTAP_GENERATE_TRANSCRIPTOME);
                        }
                        break;
                    }

                    case SIMILARITY_SEARCH: {
                        FS_dprint("STATE - SIMILARITY SEARCH");
                        // Spawn sim search object
                        std::unique_ptr<SimilaritySearch> sim_search(new SimilaritySearch(
                                databases,
                                input_path,
                                entap_data_ptrs
                        ));
                        sim_search->execute();
                        pQUERY_DATA->set_is_success_sim_search(true);
                        break;
                    }

                    case GENE_ONTOLOGY: {
                        FS_dprint("STATE - GENE ONTOLOGY");
                        std::unique_ptr<Ontology> ontology(new Ontology(
                                input_path,
                                entap_data_ptrs
                        ));
                        ontology->execute();
                        break;
                    }
                    default:
                        executeStates = EXIT;
                        break;
                }   // END SWITCH

                // Transition "states" to the next
                verify_state(state_queue, state_flag);

            } // END WHILE

            // *************************** Exit Stuff ********************** //
            pQUERY_DATA->final_statistics(final_out_dir);
           // pFileSystem->directory_iterate(FileSystem::FILE_ITER_DELETE_EMPTY, mOutpath);   // Delete empty files
            SAFE_DELETE(pQUERY_DATA);
            SAFE_DELETE(pGraphing_Manager);
            SAFE_DELETE(pEntap_Database);
        } catch (const ExceptionHandler &e) {
            SAFE_DELETE(pQUERY_DATA);
            SAFE_DELETE(pGraphing_Manager);
            SAFE_DELETE(pEntap_Database);
            exit_error(executeStates);
            throw e;
        }
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
    void verify_state(std::queue<char> &queue, bool test) {
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
     * Function bool valid_state(ExecuteStates state)
     *
     * Description          - Ensures computed state from verify_state() func
     *                        is valid and within the range
     *
     * Notes                - None
     *
     * @param state         - State to verify
     * @return              - None
     * ======================================================================
     */
    bool valid_state(ExecuteStates current_state) {
        return (current_state > INIT && current_state <= EXIT);
    }


    /**
     * ======================================================================
     * Function void exit_error(ExecuteStates exiting_state)
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
    void exit_error(ExecuteStates exiting_state) {
        std::stringstream ss;

        ss << "----------------------------------------------------\n";
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
            case COPY_FINAL_TRANSCRIPTOME:
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
                   "\t\tbe in the directory for whatever stage you failed\n"
                   "\t3. Check the debug.txt file that is printed after execution\n"
                   "\t4. Ensure your paths/inputs are correct in entap_config.ini\n"
                   "\t\tand log_file.txt (this will show your inputs)\n";
        ss << "----------------------------------------------------";
        std::cerr<<ss.str()<<std::endl;
    }
}
