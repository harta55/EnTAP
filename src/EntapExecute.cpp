//
// Created by harta on 3/4/17.
//

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
#include "Ontology.h"
#include "GraphingManager.h"
#include "DatabaseHelper.h"
#include <string>
#include <boost/regex.hpp>
#include <queue>
#include <iomanip>
//**************************************************************


namespace entapExecute {


    //*********************** Globals *****************************

    ExecuteStates           state;
    std::string             _outpath;
    std::string             _entap_outpath;
    short                   _ontology_flag;
    bool                    _isProtein;
    bool                    _EXPRESSION_SUCCESS;     // True if this stage was ran/success
    bool                    _FRAME_SELETION_SUCCESS;
    bool                    _SIM_SEARCH_SUCCESS;
    bool                    _ONTOLOGY_SUCCESS;

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
        print_debug("enTAP Executing...");

        std::map<std::string, QuerySequence>    SEQUENCE_MAP;    // Maintained/updated throughout
        std::vector<std::string>                databases;       // NCBI+UNIPROT+Other
        std::vector<std::string>                other_databases; // -d Command databases
        std::pair<std::string,std::string>      diamond_pair;    // best_hits.fa,no_hits.fa
        std::string                             input_path;      // FASTA changes depending on execution
        std::string                             no_database_hits;// No DIAMOND
        std::queue<char>                        state_queue;
        bool                                    state_flag;
        bool                                    is_complete;    // All input sequences are complete genes
        bool                                    blastp;         // false for blastx, true for blastp
        int                                     threads;

        state = INIT;
        state_flag = false;
        blastp = false;
        _EXPRESSION_SUCCESS = false;
        _FRAME_SELETION_SUCCESS = false;
        _SIM_SEARCH_SUCCESS = false;
        _ONTOLOGY_SUCCESS = false;

        input_path = user_input[ENTAP_CONFIG::INPUT_FLAG_TRANSCRIPTOME].as<std::string>();
        threads = get_supported_threads(user_input);
        _isProtein = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_RUNPROTEIN);
        is_complete = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_COMPLETE);
        diamond_pair = std::make_pair(input_path,"");

        boostFS::path working_dir(boostFS::current_path());
        _outpath = (working_dir / boostFS::path(user_input["tag"].as<std::string>())).string();
        _entap_outpath = (boostFS::path(_outpath) / boostFS::path(ENTAP_OUTPUT)).string();
        boostFS::create_directories(_entap_outpath);
        boostFS::create_directories(_outpath);

        // init_databases
        if (user_input.count("database")) {
            other_databases = user_input["database"].as<std::vector<std::string>>();
        } else other_databases.push_back(ENTAP_CONFIG::NCBI_NULL);

        // init_state_control
        if (user_input.count("state")) {
            std::string user_state_str = user_input["state"].as<std::string>();
            for (char c : user_state_str) {
                state_queue.push(c);
            }
        }

        try {
            databases = verify_databases(user_input["uniprot"].as<std::vector<std::string>>(),
                                         user_input["ncbi"].as<std::vector<std::string>>(),
                                         other_databases, "");
            verify_state(state_queue, state_flag);
            SEQUENCE_MAP = init_sequence_map(input_path, is_complete);
            GraphingManager graphingManager = GraphingManager(GRAPHING_EXE);
            FrameSelection genemark = FrameSelection(input_path,_outpath, user_input, &graphingManager);
            ExpressionAnalysis rsem = ExpressionAnalysis(input_path, threads,
                                                         _outpath, user_input, &graphingManager);
            SimilaritySearch diamond = SimilaritySearch(databases, input_path, threads,
                                                        _outpath,user_input, &graphingManager);
            Ontology ontology = Ontology(threads,_outpath,input_path, user_input, &graphingManager);

            while (state != EXIT) {
                switch (state) {
                    case FRAME_SELECTION:
                        print_debug("STATE - FRAME SELECTION");
                        if (_isProtein) {
                            print_debug("Protein sequences input, skipping frame selection");
                            flag_transcripts(FRAME_SELECTION, SEQUENCE_MAP);
                        } else {
                            input_path = genemark.execute(input_path,SEQUENCE_MAP);
                            _FRAME_SELETION_SUCCESS = true;
                        }
                        blastp = true;
                        break;
                    case RSEM:
                        print_debug("STATE - EXPRESSION");
                        if (!user_input.count(ENTAP_CONFIG::INPUT_FLAG_ALIGN)) {
                            print_debug("No alignment file specified, skipping expression analysis");
                            flag_transcripts(RSEM, SEQUENCE_MAP);
                        } else {
                            input_path = rsem.execute(input_path,SEQUENCE_MAP);
                            _EXPRESSION_SUCCESS = true;
                        }
                        break;
                    case FILTER:
                        input_path = filter_transcriptome(input_path);
                        break;
                    case DIAMOND_RUN:
                        print_debug("STATE - SIM SEARCH RUN");
                        diamond.execute(input_path, blastp);
                        _SIM_SEARCH_SUCCESS = true;
                        break;
                    case DIAMOND_PARSE:
                        print_debug("STATE - SIM SEARCH PARSE");
                        diamond_pair = diamond.parse_files(input_path,SEQUENCE_MAP);
                        input_path = diamond_pair.first;
                        no_database_hits = diamond_pair.second;
                        break;
                    case GENE_ONTOLOGY:
                        print_debug("STATE - GENE ONTOLOGY");
                        ontology.execute(SEQUENCE_MAP,input_path,no_database_hits);
                        _ONTOLOGY_SUCCESS = true;
                        break;
                    default:
                        state = EXIT;
                        break;
                }
                verify_state(state_queue, state_flag);
            }
            final_statistics(SEQUENCE_MAP);

        } catch (const ExceptionHandler &e) {throw e;}
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
     * Function std::map<std::string, QuerySequence> init_sequence_map
     *                                  (std::string     &input_file,
                                         bool            is_complete)
     *
     * Description          - Parses input transcriptome and converts to map of
     *                        each query sequence
     *                      - This map is passed throughout EnTAP execution and
     *                        updated
     *
     * Notes                - None
     *
     * @param input_file    - Path to input transcriptome
     * @param is_complete   - Flag from user if the entire transcriptome is a
     *                        complete gene
     * @return              - None
     *
     * =====================================================================
     */
    std::map<std::string, QuerySequence> init_sequence_map(std::string &input_file,
                                                           bool is_complete) {
        print_debug("Processing transcriptome...");

        std::stringstream                        out_msg;
        std::map<std::string, QuerySequence>     seq_map;
        std::string                              out_name;
        std::string                              out_new_path;
        std::string                              line;
        std::string                              sequence;
        std::string                              seq_id;
        std::string                              longest_seq;
        std::string                              shortest_seq;
        unsigned long                            count_seqs=0;
        unsigned long                            total_len=0;
        unsigned long                            shortest_len=10000;
        unsigned long                            longest_len=0;
        double                                   avg_len;
        std::vector<unsigned long>               sequence_lengths;
        std::pair<unsigned long, unsigned long>  n_vals;

        if (!file_exists(input_file)) {
            throw ExceptionHandler("Input file not found at: " + input_file,ENTAP_ERR::E_INPUT_PARSE);
        }

        boostFS::path path(input_file);
        out_name = path.filename().string();
        out_new_path = (boostFS::path(_outpath)/ boostFS::path(ENTAP_OUTPUT)).string() + out_name;
        boostFS::remove(out_new_path);
        std::ifstream in_file(input_file);
        std::ofstream out_file(out_new_path,std::ios::out | std::ios::app);

        while (true) {
            std::getline(in_file, line);
            if (line.empty() && !in_file.eof()) continue;
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            out_file << line << std::endl;
            if (line.find(">") == 0 || in_file.eof()) {
                if (!seq_id.empty()) {
                    QuerySequence query_seq = QuerySequence(_isProtein,sequence);
                    if (is_complete) query_seq.setFrame("Complete");
                    query_seq.setQseqid(seq_id);
                    seq_map.emplace(seq_id, query_seq);
                    count_seqs++;
                    unsigned long len = query_seq.getSeq_length();
                    total_len += len;
                    if (len > longest_len) {
                        longest_len = len;longest_seq = seq_id;
                    }
                    if (len < shortest_len) {
                        shortest_len = len;shortest_seq = seq_id;
                    }
                    sequence_lengths.push_back(len);
                }
                if (in_file.eof()) break;
                seq_id = line.substr(line.find(">")+1);
                sequence = line + "\n";
            } else {
                sequence += line + "\n";
            }
        }
        in_file.close(); out_file.close();
        avg_len = total_len / count_seqs;
        // first - n50, second - n90
        n_vals = calculate_N_vals(sequence_lengths, total_len);

        out_msg<<std::fixed<<std::setprecision(2);
        out_msg << ENTAP_STATS::SOFTWARE_BREAK
                << "Transcriptome Statistics\n"
                << ENTAP_STATS::SOFTWARE_BREAK<<
                "Total sequences: "                            << count_seqs    <<
                "\nTotal length of transcriptome(bp): "        << total_len     <<
                "\nAverage sequence length(bp): "              << avg_len       <<
                "\nn50: "                                      << n_vals.first  <<
                "\nn90: "                                      << n_vals.second <<
                "\nLongest sequence(bp): " << longest_len << " ("<<longest_seq<<")"<<
                "\nShortest sequence(bp): "<< shortest_len<<" ("<<shortest_seq<<")";
        if (is_complete)out_msg<<"\nAll sequences ("<<count_seqs<<") were flagged as complete genes";
        std::string msg = out_msg.str();
        print_statistics(msg);
        print_debug("Success!");
        input_file = out_new_path;
        return seq_map;
    }


    /**
     * Description - This function calculates n50 and n90 values with sequence
     *               length (nucleotide) information from a transcriptome
     *
     * @param seq_lengths - Vector of all (nucl) sequence lengths in transcriptome.
     *                    These will be sorted.
     * @param total_len   - Sum of all nucleotide lengths in transcriptome
     * @return            - Pair of <n50,n90>
     */
    std::pair<unsigned long, unsigned long> calculate_N_vals
            (std::vector<unsigned long> &seq_lengths, unsigned long total_len) {
        std::sort(seq_lengths.begin(),seq_lengths.end());
        unsigned long temp_len=0, n_50=0,n_90=0;
        double fifty_len = total_len * 0.5;
        double ninety_len = total_len * 0.9;
        for (unsigned long val : seq_lengths) {
            temp_len += val;
            if (temp_len > fifty_len && n_50 == 0) n_50 = val;
            if (temp_len > ninety_len) {
                n_90 = val;
                break;
            }
        }
        return std::pair<unsigned long, unsigned long> (n_50,n_90);
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
            state = static_cast<ExecuteStates>(state + 1);
            if (!valid_state(state)) state = EXIT;
            return;
        }
        char first = queue.front();
        if (first == 'x') {
            queue.pop();
            test = false;
            if (queue.empty()) {
                state = EXIT;
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
                state = static_cast<ExecuteStates>(state + 1);
                if (!valid_state(state)) state = EXIT;
                return;
            }
            verify_state(queue, test);
        } else {
            // some number (assuming it was parsed before)
            int i = first - '0';
            if (state < i && test) {
                state = static_cast<ExecuteStates>(state + 1);
                return;
            } else if (state < i && !test) {
                state = static_cast<ExecuteStates>(i);
                return;
            }
            if (state == i) {
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
 * Function void flag_transcripts(ExecuteStates state,
 *                              std::map<std::string, QuerySequence>& map)
 *
 * Description          - Sets boolean flags if a certain stage is skipped
 *                        in pipeline
 *                      - Used for statistics
 *
 * Notes                - Will be moved to transcriptome/project object
 *
 * @param state         - Current enum state of execution
 * @param map           - Map of query sequences (keyed to query IDs)
 * @return              - None
 * ======================================================================
 */
    void flag_transcripts(ExecuteStates state, std::map<std::string, QuerySequence>& map) {
        for (auto &pair : map) {
            switch (state) {
                case RSEM:
                    pair.second.set_is_expression_kept(true);
                    break;
                case FRAME_SELECTION:
                    pair.second.setIs_protein(true);        // Probably already done
                    break;
                default:
                    break;
            }
        }
    }


/**
 * ======================================================================
 * Function final_statistics(std::map<std::string, QuerySequence> &SEQUENCE_MAP)
 *
 * Description          - Calculates final statistical information after
 *                        completed execution
 *                      - Compiles stats on each stage of pipeline
 *
 * Notes                - None
 *
 * @param SEQUENCE_MAP  - Map of each query sequence + data
 *
 * @return              - None
 *
 * =====================================================================
 */
    void final_statistics(std::map<std::string, QuerySequence> &SEQUENCE_MAP) {
        print_debug("Pipeline finished! Calculating final statistics...");

        std::stringstream      ss;
        unsigned int           count_total_sequences=0;
        unsigned int           count_exp_kept=0;
        unsigned int           count_exp_reject=0;
        unsigned int           count_frame_kept=0;
        unsigned int           count_frame_rejected=0;
        unsigned int           count_sim_hits=0;
        unsigned int           count_sim_no_hits=0;
        unsigned int           count_ontology=0;
        unsigned int           count_no_ontology=0;
        unsigned int           count_one_go=0;
        unsigned int           count_one_kegg=0;
        unsigned int           count_sim_only=0;
        unsigned int           count_ontology_only=0;
        unsigned int           count_TOTAL_ann=0;
        unsigned int           count_TOTAL_unann=0;
        std::string            out_unannotated_nucl_path;
        std::string            out_unannotated_prot_path;
        std::string            out_annotated_nucl_path;
        std::string            out_annotated_prot_path;
        std::string            out_msg;

        out_unannotated_nucl_path = (boostFS::path(_outpath) / boostFS::path(OUT_UNANNOTATED_NUCL)).string();
        out_unannotated_prot_path = (boostFS::path(_outpath) / boostFS::path(OUT_UNANNOTATED_PROT)).string();
        out_annotated_nucl_path   = (boostFS::path(_outpath) / boostFS::path(OUT_ANNOTATED_NUCL)).string();
        out_annotated_prot_path   = (boostFS::path(_outpath) / boostFS::path(OUT_ANNOTATED_PROT)).string();

        // Switch to entap_out dir
        boostFS::remove(out_unannotated_nucl_path);
        boostFS::remove(out_unannotated_prot_path);
        boostFS::remove(out_annotated_nucl_path);
        boostFS::remove(out_annotated_prot_path);

        std::ofstream file_unannotated_nucl(out_unannotated_nucl_path, std::ios::out | std::ios::app);
        std::ofstream file_unannotated_prot(out_unannotated_prot_path, std::ios::out | std::ios::app);
        std::ofstream file_annotated_nucl(out_annotated_nucl_path, std::ios::out | std::ios::app);
        std::ofstream file_annotated_prot(out_annotated_prot_path, std::ios::out | std::ios::app);

        for (auto &pair : SEQUENCE_MAP) {
            count_total_sequences++;
            bool is_exp_kept = pair.second.is_is_expression_kept();
            bool is_prot = pair.second.isIs_protein();
            bool is_hit = pair.second.is_is_database_hit();
            bool is_ontology = pair.second.is_is_family_assigned(); // TODO Fix for interpro
            bool is_one_go = pair.second.is_is_one_go();
            bool is_one_kegg = pair.second.is_is_one_kegg();

            is_exp_kept ? count_exp_kept++ : count_exp_reject++;
            is_prot ? count_frame_kept++ : count_frame_rejected++;
            is_hit ? count_sim_hits++ : count_sim_no_hits++;
            is_ontology ? count_ontology++ : count_no_ontology++;
            if (is_one_go) count_one_go++;
            if (is_one_kegg) count_one_kegg++;

            if (is_hit && !is_ontology) count_sim_only++;
            if (!is_hit && is_ontology) count_ontology_only++;

            if (is_hit || is_ontology) {
                // Is annotated
                count_TOTAL_ann++;
                if (!pair.second.get_sequence_n().empty())
                    file_annotated_nucl<<pair.second.get_sequence_n()<<std::endl;
                if (!pair.second.get_sequence_p().empty()) {
                    file_annotated_prot<<pair.second.get_sequence_p()<<std::endl;
                }
            } else {
                // Not annotated
                if (!pair.second.get_sequence_p().empty())
                    file_unannotated_nucl<<pair.second.get_sequence_p()<<std::endl;
                if (!pair.second.get_sequence_n().empty()) {
                    file_unannotated_prot<<pair.second.get_sequence_p()<<std::endl;
                }
                count_TOTAL_unann++;
            }
        }

        file_unannotated_nucl.close();
        file_unannotated_prot.close();
        file_annotated_nucl.close();
        file_annotated_prot.close();

        ss <<
           ENTAP_STATS::SOFTWARE_BREAK          <<
           "Final Annotation Statistics\n"      <<
           ENTAP_STATS::SOFTWARE_BREAK          <<
           "Total Sequences: "                  << count_total_sequences;

        if (_EXPRESSION_SUCCESS) {
            ss <<
               "\nExpression Analysis" <<
               "\n\tKept sequences: "  << count_exp_kept    <<
               "\n\tLost sequences: "  << count_exp_reject;
        }
        if (_FRAME_SELETION_SUCCESS) {
            ss <<
               "\nFrame Selection"              <<
               "\n\tTotal sequences retained: " << count_frame_kept     <<
               "\n\tTotal sequences removed: "  << count_frame_rejected;
        }
        if (_SIM_SEARCH_SUCCESS) {
            ss <<
               "\nSimilarity Search"                               <<
               "\n\tTotal unique sequences with an alignment: "    << count_sim_hits <<
               "\n\tTotal unique sequences without an alignment: " << count_sim_no_hits;
        }
        if (_ONTOLOGY_SUCCESS) {
            switch (_ontology_flag) {
                case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
                    ss <<
                       "\nGene Families"        <<
                       "\n\tTotal unique sequences with family assignment: "    << count_ontology   <<
                       "\n\tTotal unique sequences without family assignment: " << count_no_ontology<<
                       "\n\tTotal unique sequences with at least one GO term: " << count_one_go     <<
                       "\n\tTotal unique sequences with at least one pathway (KEGG) assignment: "   << count_one_kegg;
                    break;
                case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
                    break;
                default:
                    break;
            }
        }
        ss <<
           "\nTotals"   <<
           "\n\tTotal unique sequences annotated (similarity search alignments only): "      << count_sim_only      <<
           "\n\tTotal unique sequences annotated (gene family assignment only): "            << count_ontology_only <<
           "\n\tTotal unique sequences annotated (gene family and/or similarity search): "   << count_TOTAL_ann     <<
           "\n\tTotal unique sequences unannotated (gene family and/or similarity search): " << count_TOTAL_unann;

        out_msg = ss.str();
        print_statistics(out_msg);
    }
}
