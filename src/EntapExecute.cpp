//
// Created by harta on 3/4/17.
//

#include <boost/serialization/unordered_map.hpp>
#include <iostream>
#include <boost/archive/binary_iarchive.hpp>
#include <fstream>
#include <map>
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapConsts.h"
#include "EntapInit.h"
#include <thread>
#include <list>
#include "csv.h"
#include "QuerySequence.h"
#include "FrameSelection.h"
#include "ExpressionAnalysis.h"
#include "SimilaritySearch.h"
#include "Ontology.h"
#include <string>
#include <boost/regex.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <queue>
#include <iomanip>

namespace boostFS = boost::filesystem;

namespace entapExecute {
    ExecuteStates           state;
    std::string             _frame_selection_exe;
    std::string             _expression_exe;
    std::string             _diamond_exe;
    std::string             _outpath;
    std::string             _entap_outpath;
    std::string             _ontology_exe;
    std::string             _log_path;
    short                   _ontology_flag;
    bool                    _isProtein;

    void execute_main(boost::program_options::variables_map &user_input, std::string exe_path,
                      std::unordered_map<std::string, std::string> &config_map) {
        entapInit::print_msg("enTAP Executing...");

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

        input_path = user_input["input"].as<std::string>();
        threads = entapInit::get_supported_threads(user_input);
        _isProtein = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_RUNPROTEIN);
        is_complete = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_COMPLETE);
        diamond_pair = std::make_pair(input_path,"");

        boostFS::path working_dir(boostFS::current_path());
        _outpath = working_dir.string() + "/" + user_input["tag"].as<std::string>() + "/";
        _entap_outpath = _outpath + ENTAP_EXECUTE::ENTAP_OUTPUT;
        boostFS::create_directories(_entap_outpath);
        _log_path = _outpath + ENTAP_CONFIG::LOG_FILENAME;

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
                                         other_databases, exe_path, config_map);
            init_exe_paths(config_map, exe_path);
            verify_state(state_queue, state_flag);
            SEQUENCE_MAP = init_sequence_map(input_path, is_complete);

            FrameSelection genemark = FrameSelection(input_path, _frame_selection_exe,
                                                     _outpath, user_input);
            ExpressionAnalysis rsem = ExpressionAnalysis(input_path, threads, _expression_exe,
                                                         _outpath, user_input);
            SimilaritySearch diamond = SimilaritySearch(databases, input_path, threads,
                                                        _diamond_exe, _outpath, exe_path,user_input);
            Ontology ontology = Ontology(threads,_ontology_exe,_outpath,exe_path,input_path,
                                         user_input);

            while (state != EXIT) {
                switch (state) {
                    case FRAME_SELECTION:
                        entapInit::print_msg("STATE - FRAME SELECTION");
                        if (_isProtein) {
                            entapInit::print_msg("Protein sequences input, skipping frame selection");
                        } else input_path = genemark.execute(input_path,SEQUENCE_MAP);
                        blastp = true;
                        break;
                    case RSEM:
                        entapInit::print_msg("STATE - EXPRESSION");
                        if (!user_input.count(ENTAP_CONFIG::INPUT_FLAG_ALIGN)) {
                            entapInit::print_msg("No alignment file specified, skipping expression analysis");
                        } else input_path = rsem.execute(input_path,SEQUENCE_MAP);
                        break;
                    case FILTER:
                        input_path = filter_transcriptome(input_path);
                        break;
                    case DIAMOND_RUN:
                        entapInit::print_msg("STATE - SIM SEARCH RUN");
                        diamond.execute(input_path, blastp);
                        break;
                    case DIAMOND_PARSE:
                        entapInit::print_msg("STATE - SIM SEARCH PARSE");
                        diamond_pair = diamond.parse_files(input_path,SEQUENCE_MAP);
                        input_path = diamond_pair.first;
                        no_database_hits = diamond_pair.second;
                        break;
                    case GENE_ONTOLOGY:
                        entapInit::print_msg("STATE - GENE ONTOLOGY");
                        ontology.execute(SEQUENCE_MAP,input_path,no_database_hits);
                        break;
                    default:
                        state = EXIT;
                        break;
                }
                verify_state(state_queue, state_flag);
            }
        } catch (ExceptionHandler &e) {throw e; }
    }


    std::vector<std::string> verify_databases(std::vector<std::string> uniprot, std::vector<std::string> ncbi,
                                            std::vector<std::string> database, std::string &exe,
                                            std::unordered_map<std::string, std::string> &config) {
        entapInit::print_msg("Verifying databases...");
        // return file paths
        // config file paths already exist (checked in main)
        std::vector<std::string>        file_paths;
        std::string                     path;
        std::string                     config_path;

        entapInit::print_msg("Verifying uniprot databases...");
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
                        entapInit::print_msg("Config file database found, using this path at: " +
                                             config_path);
                        path = config_path;
                    } else {
                        path = exe + ENTAP_CONFIG::UNIPROT_INDEX_PATH + u_flag + ".dmnd";
                    }
                    if (!entapInit::file_exists(path))
                        throw ExceptionHandler("Database located at: " + path + " not found", ENTAP_ERR::E_INPUT_PARSE);
                    file_paths.push_back(path);
                } else {
                    entapInit::print_msg("No/null Uniprot databases detected");
                    break;
                }
            }
        }
        entapInit::print_msg("Complete");
        entapInit::print_msg("Verifying NCBI databases...");
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
                        entapInit::print_msg("Config file database found, using this path at: " +
                                             config_path);
                        path = config_path;
                    } else {
                        path = exe + ENTAP_CONFIG::NCBI_INDEX_PATH + u_flag + ".dmnd";
                    }
                    if (!entapInit::file_exists(path))
                        throw ExceptionHandler("Database located at: " + path + " not found", ENTAP_ERR::E_INPUT_PARSE);
                    file_paths.push_back(path);
                } else {
                    entapInit::print_msg("No/null NCBI databases detected");
                    break;
                }
            }
        }
        entapInit::print_msg("Complete");
        entapInit::print_msg("Verifying other databases...");
        if (database.size() > 0) {
            for (auto const &data_path:database) {
                if (data_path.compare(ENTAP_CONFIG::NCBI_NULL) == 0) continue;
                if (!entapInit::file_exists(data_path)) {
                    throw ExceptionHandler("Database located at: " + data_path + " not found",
                                           ENTAP_ERR::E_INPUT_PARSE);
                }
                boostFS::path bpath(data_path);
                std::string ext = bpath.extension().string();
                if (ext.compare(".dmnd") == 0) {
                    entapInit::print_msg("User has input a diamond indexed database at: " +
                                         data_path);
                    file_paths.push_back(data_path);
                    continue;
                } else {
                    //todo fix
                    entapInit::print_msg("User has input a database at: " + data_path);
                    std::string test_path = exe + ENTAP_CONFIG::BIN_PATH + data_path + ".dmnd";
                    entapInit::print_msg("Checking if indexed file exists at: " + test_path);
                    if (!entapInit::file_exists(test_path)) {
                        throw ExceptionHandler("Database located at: " + data_path + " not found",
                                               ENTAP_ERR::E_INPUT_PARSE);
                    } else {
                        file_paths.push_back(test_path);
                    }
                }
            }
        }
        entapInit::print_msg("Verification complete!");
        if (file_paths.size() > 0) {
            std::string database_final = "\n\nDatabases selected:\n";
            for (std::string base: file_paths) {
                database_final += base + "\n";
            }
            entapInit::print_msg(database_final);
        } else {
            entapInit::print_msg("No databases selected, some funcionality "
                                         "may not be able to run");
        }
        return file_paths;
    }


    /**
     * Description          - Merely selects transcriptome that will continue in pipeline
     *                        and copies it to entap_out directory
     *
     * @param input_path    - Input transcriptome (expression and/or frame selected)
     * @return              - Copied transcriptome
     */
    std::string filter_transcriptome(std::string &input_path) {
        entapInit::print_msg("Beginning to copy final transcriptome to be used...");
        boostFS::path file_name(input_path);
        file_name.filename();
        std::string out_path = _entap_outpath +
                               file_name.stem().stem().string() + "_final.fasta";
        boostFS::copy_file(input_path,out_path,boostFS::copy_option::overwrite_if_exists);
        return out_path;
    }


    std::map<std::string, QuerySequence> init_sequence_map(std::string &input_file,
                                                           bool is_complete) {
        std::stringstream out_msg;out_msg<<std::fixed<<std::setprecision(2);
        entapInit::print_msg("Processing transcriptome...");

        std::map<std::string, QuerySequence>            seq_map;



        if (!entapInit::file_exists(input_file)) {
            throw ExceptionHandler("Input file not found at: " +
                input_file,ENTAP_ERR::E_INPUT_PARSE);
        }
        out_msg << ENTAP_STATS::SOFTWARE_BREAK
                << "Transcriptome Statistics\n"
                << ENTAP_STATS::SOFTWARE_BREAK;
        boostFS::path path(input_file);
        std::string out_name = path.filename().string();
        std::string out_new_path = _outpath + ENTAP_EXECUTE::ENTAP_OUTPUT + out_name;
        boostFS::remove(out_new_path);
        std::ifstream in_file(input_file);
        std::ofstream out_file(out_new_path,std::ios::out | std::ios::app);
        std::string line, sequence, seq_id, longest_seq, shortest_seq;
        unsigned long count_seqs=0, total_len=0,shortest_len = 10000, longest_len = 0;
        std::vector<unsigned long> sequence_lengths;
        while (true) {
            std::getline(in_file, line);
            if (line.empty() && !in_file.eof()) continue;
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            out_file << line << std::endl;
            if (line.find(">") == 0 || in_file.eof()) {
                if (!seq_id.empty()) {
                    QuerySequence query_seq = QuerySequence(_isProtein,sequence);
                    if (is_complete) query_seq.setFrame(ENTAP_EXECUTE::FRAME_SELECTION_COMPLETE_FLAG);
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
        double avg_len = total_len / count_seqs;
        // first - n50, second - n90
        std::pair<unsigned long, unsigned long> n_vals =
                calculate_N_vals(sequence_lengths, total_len);
        out_msg <<
                "Total sequences: "                            << count_seqs    <<
                "\nTotal length of transcriptome(bp): "        << total_len     <<
                "\nAverage sequence length(bp): "              << avg_len       <<
                "\nn50: "                                      << n_vals.first  <<
                "\nn90: "                                      << n_vals.second <<
                "\nLongest sequence(bp): " << longest_len << " ("<<longest_seq<<")"<<
                "\nShortest sequence(bp): "<< shortest_len<<" ("<<shortest_seq<<")";
        if (is_complete)out_msg<<"\nAll sequences ("<<count_seqs<<") were flagged as complete genes";
        std::string msg = out_msg.str();
        print_statistics(msg,_outpath);
        entapInit::print_msg("Success!");
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
     * Description      - Handles log file information output to .txt file.
     *
     * @param msg       - String message
     * @param out_path  - enTAP outfiles directory
     */
    void print_statistics(std::string &msg, std::string &out_path) {
        std::string file_path = out_path + ENTAP_CONFIG::LOG_FILENAME;
        std::ofstream log_file(file_path, std::ios::out | std::ios::app);
        log_file << msg << std::endl;
        log_file.close();
    }


    /**
     * Description      - Finds exe paths for each piece of pipeline. WARNING: does not
     *                    check for validity of exe paths (as some parts may not want to be
     *                    ran)
     *
     * @param map       - Map of entap_config file
     * @param exe       - enTAP execution directory
     * @return          - DIAMOND .exe path ran by enTAP::Init
     */
    std::string init_exe_paths(std::unordered_map<std::string, std::string> &map, std::string &exe) {
        entapInit::print_msg("Verifying execution paths. Note they are not checked for validity...");
        std::string temp_rsem       = map[ENTAP_CONFIG::KEY_RSEM_EXE];
        std::string temp_diamond    = map[ENTAP_CONFIG::KEY_DIAMOND_EXE];
        std::string temp_genemark   = map[ENTAP_CONFIG::KEY_GENEMARK_EXE];
        std::string temp_eggnog     = map[ENTAP_CONFIG::KEY_EGGNOG_EXE];
        std::string temp_interpro   = map[ENTAP_CONFIG::KEY_INTERPRO_EXE];

        if (_ontology_flag>1 || _ontology_flag<0) {
            throw ExceptionHandler("Annotation flag must be either 0(eggnog) or"
                                           "1(interproscan)", ENTAP_ERR::E_CONFIG_PARSE);}
        if (temp_rsem.empty()) {
            temp_rsem = exe + ENTAP_EXECUTE::RSEM_EXE_PATH;
            entapInit::print_msg("RSEM config path empty, setting to default: " + temp_rsem);
        } else entapInit::print_msg("RSEM path set to: " + temp_rsem);

        if (temp_diamond.empty()) {
            temp_diamond = exe + ENTAP_CONFIG::DIAMOND_PATH_EXE;
            entapInit::print_msg("DIAMOND config path empty, setting to default: "+temp_diamond);
        } else entapInit::print_msg("DIAMOND path set to: " + temp_diamond);

        if (temp_genemark.empty()) {
            temp_genemark = exe + ENTAP_EXECUTE::GENEMARK_EXE_PATH;
            entapInit::print_msg("GenemarkS-T config path empty, setting to default: "+temp_genemark);
        } else entapInit::print_msg("GenemarkS-T path set to: " + temp_genemark);

        if (temp_eggnog.empty()) {
            temp_eggnog = exe + ENTAP_EXECUTE::EGGNOG_EMAPPER_EXE;
            entapInit::print_msg("Eggnog config path empty, setting to default: " + temp_eggnog);
        } else entapInit::print_msg("Eggnog path set to: " + temp_eggnog);

        if (temp_interpro.empty()) {
            temp_interpro = exe + ENTAP_EXECUTE::INTERPRO_EXE;
            entapInit::print_msg("Interpro config path empty, setting to default: "+temp_interpro);
        } else entapInit::print_msg("Interpro path set to: " + temp_interpro);

        _diamond_exe = temp_diamond;
        _frame_selection_exe = temp_genemark;
        _expression_exe = temp_rsem;
        _ontology_flag==0 ? _ontology_exe = temp_eggnog : _ontology_exe = temp_interpro;
        entapInit::print_msg("Success! All exe paths set");
        return _diamond_exe;        // for config run
    }

    //only assuming between 0-9 NO 2 DIGIT STATES
    void verify_state(std::queue<char> &queue, bool &test) {
        entapInit::print_msg("verifying state...");
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
        entapInit::print_msg("Success!");
    }

    bool valid_state(ExecuteStates s) {
        return (s >= RSEM && s <= EXIT);
    }
}
