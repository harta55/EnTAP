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
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <csv.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/regex.hpp>
#include <iomanip>
#include "SimilaritySearch.h"
#include "ExceptionHandler.h"
#include "EntapExecute.h"
#include "EntapGlobals.h"
#include "common.h"
//**************************************************************


/**
 * ======================================================================
 * Function SimilaritySearch(std::vector<std::string> &databases, std::string input,
                           int threads, std::string out, boost::program_options::variables_map &user_flags,
                           GraphingManager *graphingManager)
 *
 * Description          - SimilaritySearch object constructor
 *                      - Responsible for initiating SimSearch member variables,
 *                        parsing user input for relevant information used within
 *                        this module
 *
 * Notes                - None
 *
 * @param databases     - List of user selected databases (parsed in  execute namespace)
 * @param input         - Path to user transcriptome (final version post filtering)
 * @param threads       - Thread count
 * @param out           - EnTAP out directory
 * @param user_flags    - Boost parsed user input
 * @param graphingManager- Pointer to graphing manager
 *
 * @return              - SimilaritySearch instance
 * ======================================================================
 */
SimilaritySearch::SimilaritySearch(std::vector<std::string> &databases, std::string input,
                           int threads, std::string out, boost::program_options::variables_map &user_flags,
                           GraphingManager *graphingManager, QueryData *queryData) {
    print_debug("Spawn object - SimilaritySearch");
    _pQUERY_DATA    = queryData;
    _database_paths = databases;
    _input_path     = input;
    _threads        = threads;
    _diamond_exe    = DIAMOND_EXE;      // Set to extern set previously
    _outpath        = out;
    _input_species  = "";

    // Species already checked for validity in Init
    if (user_flags.count(ENTAP_CONFIG::INPUT_FLAG_SPECIES)) {
        _input_species = user_flags[ENTAP_CONFIG::INPUT_FLAG_SPECIES].as<std::string>();
        std::transform(_input_species.begin(), _input_species.end(), _input_species.begin(), ::tolower);
        std::replace(_input_species.begin(), _input_species.end(), '_',' ');
    }
    _qcoverage = user_flags[ENTAP_CONFIG::INPUT_FLAG_QCOVERAGE].as<fp32>();
    _tcoverage = user_flags[ENTAP_CONFIG::INPUT_FLAG_TCOVERAGE].as<fp32>();
    _overwrite = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _e_val     = user_flags[ENTAP_CONFIG::INPUT_FLAG_E_VAL].as<fp32>();

    // Format contaminants for use in database
    std::vector<std::string> contaminants;
    if (user_flags.count(ENTAP_CONFIG::INPUT_FLAG_CONTAM)) {
        contaminants = user_flags[ENTAP_CONFIG::INPUT_FLAG_CONTAM].as<std::vector<std::string>>();
    }
    for (int ind = 0; ind < contaminants.size(); ind++) {
        if (contaminants[ind].empty()) continue;
        std::string &str = contaminants[ind];
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        std::replace(str.begin(), str.end(), '_',' ');
    }
    _contaminants    = contaminants;
    _software_flag   = DIAMOND_FLAG;         // Default DIAMOND software
    _pGraphingManager= graphingManager;

    // Set sim search paths/directories
    _sim_search_dir  = PATHS(out, SIM_SEARCH_DIR);
    _processed_path  = PATHS(_sim_search_dir, PROCESSED_DIR);
    _results_path    = PATHS(_sim_search_dir, RESULTS_DIR);
    _figure_path     = PATHS(_sim_search_dir, FIGURE_DIR);
}


/**
 * ======================================================================
 * Function SimilaritySearch()
 *
 * Description          - Empty constructor
 *
 * Notes                - None
 *
 * @return              - SimilaritySearch instance
 * ======================================================================
 */
SimilaritySearch::SimilaritySearch() {}


/**
 * ======================================================================
 * Function std::vector<std::string> SimilaritySearch::execute(std::string updated_input,
 *                                                             bool blast)
 *
 * Description          - Responsible for executing the selected similarity
 *                        searching software
 *                      - Returns output files of sim search
 *                      - Throws instance of ExceptionHandler on failed execution
 *
 * Notes                - None
 *
 * @param updated_input - User transcriptome, if changed
 * @param blast         - Blastx/blastp (true for blastp)
 *
 * @return              - Vector of output files
 * ======================================================================
 */
std::vector<std::string> SimilaritySearch::execute(std::string updated_input,bool blast) {
    this->_input_path = updated_input;
    this->_blastp = blast;
    _blastp ? _blast_type = "blastp" : _blast_type = "blastx";
    try {
        switch (_software_flag) {
            case 0:
                return diamond();
            default:
                return diamond();
        }
    } catch (ExceptionHandler &e) {throw e;}
}


/**
 * ======================================================================
 * Function std::pair<std::string,std::string> SimilaritySearch::parse_files(std::string new_input,
                                   std::map<std::string, QuerySequence>& MAP)
 *
 * Description          - Responsible for best hit selection of sequences hit
 *                        against the diamond databases
 *                      - Throws fatal ExceptionHandler instance
 *
 * Notes                - Entered as SIM_SEARCH_PARSE state from Execute namespace
 *
 * @param new_input     - User transcriptome, if changed
 * @param MAP           - Master data structure of transcriptome information
 *
 * @return              - Pair of output files (hits and no hits)
 * ======================================================================
 */
std::pair<std::string,std::string> SimilaritySearch::parse_files(std::string new_input) {
    _input_path = new_input;
    try {
        switch (_software_flag) {
            case DIAMOND_FLAG:
                return diamond_parse(_contaminants);
            default:
                return diamond_parse(_contaminants);
        }
    } catch (ExceptionHandler &e) {throw e;}
}


/**
 * ======================================================================
 * Function std::vector<std::string> SimilaritySearch::diamond()
 *
 * Description          - Responsible for executing simliarity search through
 *                        pstreams library
 *                      - Returns vector of output files from sim search
 *                      - Checks whether DIAMOND has been ran previously
 *
 * Notes                - None
 *
 *
 * @return              - Output files from similarity searching
 * ======================================================================
 */
std::vector<std::string> SimilaritySearch::diamond() {
    print_debug("Beginning to execute DIAMOND...");

    std::vector<std::string>    out_paths;
    boostFS::path               transc_name;
    std::string                 filename;
    std::string                 out_path;
    std::string                 std_out;

    if (!file_exists(_input_path)) {
        throw ExceptionHandler("Transcriptome file not found",ENTAP_ERR::E_RUN_SIM_SEARCH_RUN);
    }
    transc_name = _input_path;
    transc_name = transc_name.stem();
    if (transc_name.has_stem()) transc_name = transc_name.stem(); //.fasta.faa
    if (_overwrite) {
        boostFS::remove_all(_sim_search_dir);
    }
    boostFS::create_directories(_sim_search_dir);
    // database verification already ran, don't need to verify each path
    try {
        // assume all paths should be .dmnd
        for (std::string data_path : _database_paths) {
            print_debug("Searching against database located at: " + data_path + "...");
            boostFS::path database_name(data_path);
            database_name = database_name.stem();
            filename = _blast_type + "_" + transc_name.string() + "_" + database_name.string();
            out_path = PATHS(_sim_search_dir,filename) + ".out";
            std_out  = PATHS(_sim_search_dir,filename) + "_std";
            _file_to_database[out_path] = database_name.string();
            if (file_exists(out_path)) {
                print_debug("File found at " + out_path + " skipping execution against this database");
                out_paths.push_back(out_path);
                continue;
            }
            diamond_blast(_input_path, out_path, std_out,data_path, _threads, _blast_type);
            print_debug("Success! Results written to " + out_path);
            out_paths.push_back(out_path);
        }
    } catch (const ExceptionHandler &e) {throw e;}
    _sim_search_paths = out_paths;
    return out_paths;
}


/**
 * ======================================================================
 * Function void SimilaritySearch::diamond_blast(std::string input_file, std::string output_file, std::string std_out,
                   std::string &database,int &threads, std::string &blast)
 *
 * Description          - Responsible for execution of DIAMOND through pstreams
 *                        library
 *
 * Notes                - None
 *
 * @param input_file    - Path to input transcriptome
 * @param output_file   - Path to output file from sim search
 * @param std_out       - Std out/err path
 * @param database      - Selected database to hit against
 * @param threads       - Thread number
 * @param blast         - Blast type (blastx/blastp)
 *
 * @return              - None
 * ======================================================================
 */
void SimilaritySearch::diamond_blast(std::string input_file, std::string output_file, std::string std_out,
                   std::string &database,int &threads, std::string &blast) {

    std::string        diamond_run;

    diamond_run =
            _diamond_exe + " "
            + blast +
            " -d " + database    +
            " --query-cover "    + std::to_string(_qcoverage) +
            " --subject-cover "  + std::to_string(_tcoverage) +
            " --more-sensitive"  +
            " --top 3"           +
            " -q " + input_file  +
            " -o " + output_file +
            " -p " + std::to_string(threads) +
            " -f " + "6 qseqid sseqid pident length mismatch gapopen "
                     "qstart qend sstart send evalue bitscore qcovhsp stitle";

    print_debug("\nExecuting Diamond:\n" + diamond_run);
    if (execute_cmd(diamond_run, std_out) != 0) {
        throw ExceptionHandler("Error in DIAMOND run with database located at: " +
                               database, ENTAP_ERR::E_RUN_SIM_SEARCH_RUN);
    }
}


/**
 * ======================================================================
 * Function std::vector<std::string> SimilaritySearch::verify_diamond_files(
 *                                      std::string &outpath,
 *                                      std::string name)
 *
 * Description          - Checks whether DIAMOND has already been ran with
 *                        the same parameters to skip re-running
 *
 * Notes                - None
 *
 * @param outpath       - Path to out directory
 * @param name          - Name of input transcriptome
 *
 * @return              - None
 * ======================================================================
 */
std::vector<std::string> SimilaritySearch::verify_diamond_files(std::string &outpath, std::string name) {
    print_debug("Override unselected, checking for diamond files of selected databases...");
    std::vector<std::string> out_list;
    std::string              temp_out;
    std::string              file_name_full;

    for (std::string data_path : _database_paths) {
        // assume all paths should be .dmnd
        boostFS::path file_name(data_path);
        file_name = file_name.stem();
        file_name_full = _blast_type + "_" + name + "_" + file_name.string() + ".out";
        temp_out = PATHS(outpath, file_name_full);
        if (!file_exists(temp_out)){
            print_debug("File at: " + temp_out + " not found, running diamond");
            out_list.clear();return out_list;
        }
        out_list.push_back(temp_out);
    }
    print_debug("All diamond files found, skipping this stage of enTAP");
    _sim_search_paths = out_list;
    return out_list;
}

// input: 3 database string array of selected databases
std::pair<std::string,std::string> SimilaritySearch::diamond_parse(std::vector<std::string>& contams) {
    print_debug("Beginning to filter individual diamond_files...");

    std::unordered_map<std::string, std::string>    taxonomic_database;
    std::list<std::map<std::string,QuerySequence>>  database_maps;
    uint32                                          count_removed;
    uint32                                          count_TOTAL_hits;
    uint32                                          count_under_e;

    try {
        taxonomic_database = read_tax_map();
    } catch (ExceptionHandler &e) {throw e;}
    if (_sim_search_paths.empty()) {
        boostFS::path transc_name(_input_path); transc_name=transc_name.stem();
        if (transc_name.has_stem()) transc_name = transc_name.stem(); //.fasta.faa
        _sim_search_paths = verify_diamond_files(_sim_search_dir,transc_name.string());
    }
    if (_sim_search_paths.empty()) throw ExceptionHandler("No diamond files found", ENTAP_ERR::E_RUN_SIM_SEARCH_FILTER);

    boostFS::remove_all(_processed_path);
    boostFS::create_directories(_processed_path);
    _input_lineage = get_lineage(_input_species,taxonomic_database);

    for (std::string &data : _sim_search_paths) {
        print_debug("Diamond file located at " + data + " being filtered");
        io::CSVReader<DMND_COL_NUMBER, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(data);
        // todo have columns from input file, in_read_header for versatility
        std::string qseqid, sseqid, stitle, database_name,pident, bitscore,
                length, mismatch, gapopen, qstart, qend, sstart, send;
        fp64 evalue;
        fp64 coverage;
        count_removed    = 0;
        count_TOTAL_hits = 0;
        count_under_e    = 0;
        std::stringstream out_stream;

        if (_file_to_database.find(data) != _file_to_database.end()) {
            database_name = _file_to_database[data];
        } else database_name = boostFS::path(data).filename().stem().string();
        std::string out_base_path = (boostFS::path(_processed_path) / boostFS::path(database_name)).string();
        std::string out_unselected_tsv = (boostFS::path(out_base_path) /
                boostFS::path(SIM_SEARCH_DATABASE_UNSELECTED)).string();
        boostFS::create_directories(out_base_path);
        std::ofstream file_unselected_tsv(out_unselected_tsv,std::ios::out | std::ios::app);
        std::map<std::string, QuerySequence> database_map;
        print_header(out_unselected_tsv);

        while (in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                           qstart, qend, sstart, send, evalue, bitscore, coverage,stitle)) {
            count_TOTAL_hits++;

            QuerySequence new_query = QuerySequence();
            new_query.set_sim_search_results(data, qseqid, sseqid, pident, length, mismatch, gapopen,
                                             qstart, qend, sstart, send, stitle, bitscore, evalue,coverage);
            new_query.setIs_better_hit(true);
            std::string species = get_species(stitle);
            new_query.setSpecies(species);
            std::string lineage = get_lineage(species,taxonomic_database);
            new_query.set_lineage(lineage);
            std::pair<bool,std::string> contam_info = is_contaminant(lineage, taxonomic_database,contams);
            new_query.setContaminant(contam_info.first);
            new_query.set_contam_type(contam_info.second);
            bool informative = is_informative(stitle);
            new_query.set_is_informative(informative);
            new_query.setFrame((*_pQUERY_DATA->get_sequences_ptr()).at(qseqid).getFrame());  // May want to handle differently, SLOW
            new_query.set_tax_score(_input_lineage);
            if (evalue > _e_val) {
                count_under_e++; count_removed++;
                file_unselected_tsv << new_query.print_tsv(DEFAULT_HEADERS) <<std::endl;
                continue;
            }
            std::map<std::string, QuerySequence>::iterator it = database_map.find(qseqid);
            if (it != database_map.end()) {
                if (new_query > it->second) {
                    file_unselected_tsv << it->second.print_tsv(DEFAULT_HEADERS) << std::endl;
                    it->second = new_query;
                    count_removed++;
                } else {
                    count_removed++;
                    file_unselected_tsv << new_query.print_tsv(DEFAULT_HEADERS)<< std::endl;
                }
            } else database_map.emplace(qseqid, new_query);
        }
        print_debug("File parsed, calculating statistics and writing output...");
        // final database stats
        file_unselected_tsv.close();
        out_stream<<std::fixed<<std::setprecision(2);
        out_stream << ENTAP_STATS::SOFTWARE_BREAK
                   << "Similarity Search - Diamond - "    <<database_name<<"\n"
                   << ENTAP_STATS::SOFTWARE_BREAK <<
                   "Search results:\n"            << data <<
                   "\n\tTotal alignments: "               << count_TOTAL_hits   <<
                   "\n\tTotal unselected results: "       << count_removed      <<
                   "\n\t\tWritten to: "                   << out_unselected_tsv;
        calculate_best_stats(database_map,out_stream,out_base_path,false);
        std::string out_msg = out_stream.str() + "\n";
        print_statistics(out_msg);
        database_maps.push_back(database_map);
        print_debug("Success!");
    }
    return process_best_diamond_hit(database_maps);
}

std::pair<std::string,std::string> SimilaritySearch::calculate_best_stats (
                                             std::map<std::string, QuerySequence>&best_hits,
                                             std::stringstream &ss, std::string &base_path, bool is_final) {

    boostFS::path               base_bst;
    GraphingStruct              graphingStruct;
    std::string                 database;
    std::string                 figure_base;
    std::string                 frame;
    uint64                      count_no_hit=0;
    uint64                      count_contam=0;
    uint64                      count_filtered=0;
    uint64                      count_informative=0;
    uint64                      count_uninformative=0;
    uint32                      ct;
    fp64                        percent;
    fp64                        contam_percent;
    std::map<std::string, int>  contam_map;
    std::map<std::string, int>  species_map;
    std::map<std::string, int>  contam_species_map;
    graph_sum_t                 graphing_sum_map;

    database = boostFS::path(base_path).filename().string();
    figure_base = (boostFS::path(base_path) / FIGURE_DIR).string();
    base_bst = base_path;
    boostFS::create_directories(figure_base);
    boostFS::create_directories(base_path);     // should be created before

    std::string out_best_contams_tsv             = PATHS(base_bst ,SIM_SEARCH_DATABASE_CONTAM_TSV);
    std::string out_best_contams_fa_nucl         = PATHS(base_bst ,SIM_SEARCH_DATABASE_CONTAM_FA_NUCL);
    std::string out_best_contams_fa_prot         = PATHS(base_bst ,SIM_SEARCH_DATABASE_CONTAM_FA_PROT);

    std::string out_best_hits_tsv                = PATHS(base_bst, SIM_SEARCH_DATABASE_BEST_TSV);
    std::string out_best_hits_no_contam_tsv      = PATHS(base_bst, SIM_SEARCH_DATABASE_BEST_TSV_NO_CONTAM);
    std::string out_best_hits_fa_nucl            = PATHS(base_bst, SIM_SEARCH_DATABASE_BEST_FA_NUCL);
    std::string out_best_hits_fa_prot            = PATHS(base_bst, SIM_SEARCH_DATABASE_BEST_FA_PROT);
    std::string out_best_hits_fa_nucl_no_contam  = PATHS(base_bst, SIM_SEARCH_DATABASE_BEST_FA_NUCL_NO_CONTAM);
    std::string out_best_hits_fa_prot_no_contam  = PATHS(base_bst, SIM_SEARCH_DATABASE_BEST_FA_PROT_NO_CONTAM);

    std::string out_no_hits_fa_nucl              = PATHS(base_bst, SIM_SEARCH_DATABASE_NO_HITS_NUCL);
    std::string out_no_hits_fa_prot              = PATHS(base_bst, SIM_SEARCH_DATABASE_NO_HITS_PROT);

    std::string graph_species_txt_path           = PATHS(figure_base, GRAPH_SPECIES_BAR_TXT);
    std::string graph_species_png_path           = PATHS(figure_base, GRAPH_SPECIES_BAR_PNG);
    std::string graph_contam_txt_path            = PATHS(figure_base, GRAPH_CONTAM_BAR_TXT);
    std::string graph_contam_png_path            = PATHS(figure_base, GRAPH_CONTAM_BAR_PNG);
    std::string graph_sum_txt_path               = PATHS(figure_base, GRAPH_DATABASE_SUM_TXT);
    std::string graph_sum_png_path               = PATHS(figure_base, GRAPH_DATABASE_SUM_PNG);

    print_header(out_best_contams_tsv);
    print_header(out_best_hits_tsv);
    print_header(out_best_hits_no_contam_tsv);

    std::ofstream file_best_hits_tsv(out_best_hits_tsv,std::ios::out | std::ios::app);
    std::ofstream file_best_hits_tsv_no_contam(out_best_hits_no_contam_tsv,std::ios::out | std::ios::app);
    std::ofstream file_best_hits_fa_nucl(out_best_hits_fa_nucl,std::ios::out | std::ios::app);
    std::ofstream file_best_hits_fa_prot(out_best_hits_fa_prot,std::ios::out | std::ios::app);
    std::ofstream file_best_hits_fa_nucl_no_contam(out_best_hits_fa_nucl_no_contam,std::ios::out | std::ios::app);
    std::ofstream file_best_hits_fa_prot_no_contam(out_best_hits_fa_prot_no_contam,std::ios::out | std::ios::app);
    std::ofstream file_best_contam_tsv(out_best_contams_tsv,std::ios::out | std::ios::app);
    std::ofstream file_best_contam_fa_prot(out_best_contams_fa_prot,std::ios::out | std::ios::app);
    std::ofstream file_best_contam_fa_nucl(out_best_contams_fa_nucl,std::ios::out | std::ios::app);
    std::ofstream file_no_hits_nucl(out_no_hits_fa_nucl, std::ios::out | std::ios::app);
    std::ofstream file_no_hits_prot(out_no_hits_fa_prot, std::ios::out | std::ios::app);

    std::ofstream graph_species_file(graph_species_txt_path, std::ios::out | std::ios::app);
    std::ofstream graph_contam_file(graph_contam_txt_path, std::ios::out | std::ios::app);
    std::ofstream graph_sum_file(graph_sum_txt_path, std::ios::out | std::ios::app);

    graph_species_file << "Species\tCount"     << std::endl;
    graph_contam_file  << "Contaminant Species\tCount" << std::endl;
    graph_sum_file     << "Category\tCount"    << std::endl;

    for (auto &pair : *_pQUERY_DATA->get_sequences_ptr()) {
        std::map<std::string, QuerySequence>::iterator it = best_hits.find(pair.first);
        // Check if original sequences have hit a database
        if (it == best_hits.end()) {
            if ((pair.second.isIs_protein() && _blastp) || (!pair.second.isIs_protein() && !_blastp)) {
                // Protein/nucleotide did not hit database
                count_no_hit++;
                file_no_hits_nucl << pair.second.get_sequence_n() << std::endl;
                file_no_hits_prot << pair.second.get_sequence_p() << std::endl;
                // Graphing
                frame = pair.second.getFrame();
                if (graphing_sum_map[frame].find(NO_HIT_FLAG) != graphing_sum_map[frame].end()) {
                    graphing_sum_map[frame][NO_HIT_FLAG]++;
                } else graphing_sum_map[frame][NO_HIT_FLAG] = 1;

            }
        } else {
            // Have hit a database
            frame = pair.second.getFrame();     // Used for graphing
            file_best_hits_fa_nucl << pair.second.get_sequence_n()<<std::endl;
            file_best_hits_fa_prot << pair.second.get_sequence_p()<<std::endl;
            file_best_hits_tsv << it->second.print_tsv(DEFAULT_HEADERS) << std::endl;

            count_filtered++;
            std::string species;
            if (!it->second.get_species().empty()) {
                species = it->second.get_species();
            }
            if (it->second.isContaminant()) {
                count_contam++;
                file_best_contam_fa_nucl << pair.second.get_sequence_n()<<std::endl;
                file_best_contam_fa_prot << pair.second.get_sequence_p()<<std::endl;
                file_best_contam_tsv << it->second.print_tsv(DEFAULT_HEADERS) << std::endl;
                std::string contam = it->second.get_contam_type();
                if (contam_map.count(contam)) {
                    contam_map[contam]++;
                } else contam_map[contam] = 1;
                if (contam_species_map.count(species)) {
                    contam_species_map[species]++;
                } else contam_species_map[species] = 1;
            } else {
                file_best_hits_fa_nucl_no_contam << pair.second.get_sequence_n()<<std::endl;
                file_best_hits_fa_prot_no_contam << pair.second.get_sequence_p()<<std::endl;
                file_best_hits_tsv_no_contam << it->second.print_tsv(DEFAULT_HEADERS) << std::endl;
            }
            if (species_map.count(species)) {
                species_map[species]++;
            } else species_map[species] = 1;

            if (it->second.is_informative()) {
                count_informative++;
                // Graphing
                if (graphing_sum_map[frame].find(INFORMATIVE_FLAG) != graphing_sum_map[frame].end()) {
                    graphing_sum_map[frame][INFORMATIVE_FLAG]++;
                } else graphing_sum_map[frame][INFORMATIVE_FLAG] = 1;

            } else {
                count_uninformative++;
                if (graphing_sum_map[frame].find(UNINFORMATIVE_FLAG) != graphing_sum_map[frame].end()) {
                    graphing_sum_map[frame][UNINFORMATIVE_FLAG]++;
                } else graphing_sum_map[frame][UNINFORMATIVE_FLAG] = 1;
            }

            if (is_final) {
                pair.second.set_sim_struct(it->second.get_sim_struct());
                pair.second.set_is_database_hit(true);
            }

        }
    }
    file_best_hits_tsv.close();
    file_best_hits_tsv_no_contam.close();
    file_best_hits_fa_nucl.close();
    file_best_hits_fa_prot.close();
    file_best_hits_fa_nucl_no_contam.close();
    file_best_hits_fa_prot_no_contam.close();
    file_best_contam_tsv.close();
    file_best_contam_fa_prot.close();
    file_best_contam_fa_nucl.close();
    file_no_hits_nucl.close();
    file_no_hits_prot.close();

    std::vector<count_pair> contam_species_vect(contam_species_map.begin(),contam_species_map.end());
    std::vector<count_pair> species_vect(species_map.begin(),species_map.end());
    std::sort(contam_species_vect.begin(),contam_species_vect.end(),compair());
    std::sort(species_vect.begin(),species_vect.end(),compair());
    contam_percent = ((double)count_contam / count_filtered) * 100;

    ss <<
       "\n\tTotal unique transcripts with an alignment: "                              << count_filtered          <<
       "\n\t\tReference transcriptome sequences with an alignment (FASTA):\n\t\t\t"    << out_best_hits_fa_prot   <<
       "\n\t\tSearch results (TSV):\n\t\t\t"          << out_best_hits_tsv   <<
       "\n\tTotal unique transcripts without an alignment: "                 << count_no_hit       <<
       "\n\t\tReference transcriptome sequences without an alignment (FASTA):\n\t\t\t"    << out_no_hits_fa_prot;
    // Have frame information
    if (graphing_sum_map.size() > 1) {
        for (auto &pair : graphing_sum_map) {
            // Frame -> Map of uninform/inform/no hits
            ss << "\n\t\t" << pair.first << "(" << pair.second[NO_HIT_FLAG] << ")";
            graph_sum_file << pair.first << "\t" << NO_HIT_FLAG << "\t" << pair.second[NO_HIT_FLAG] << "\n";
        }
    }
    ss <<
       "\n\tTotal unique informative alignments: " << count_informative;
    if (graphing_sum_map.size() > 1) {
        for (auto &pair : graphing_sum_map) {
            // Frame -> Map of uninform/inform/no hits
            ss << "\n\t\t" << pair.first << "(" << pair.second[INFORMATIVE_FLAG] << ")";
            graph_sum_file << pair.first << "\t" << INFORMATIVE_FLAG << "\t" << pair.second[INFORMATIVE_FLAG] << "\n";
        }
    }
    ss <<
       "\n\tTotal unique uninformative alignments: " << count_uninformative;
    if (graphing_sum_map.size() > 1) {
        for (auto &pair : graphing_sum_map) {
            // Frame -> Map of uninform/inform/no hits
            ss << "\n\t\t" << pair.first << "(" << pair.second[UNINFORMATIVE_FLAG] << ")";
            graph_sum_file << pair.first << "\t" << UNINFORMATIVE_FLAG << "\t" << pair.second[UNINFORMATIVE_FLAG] << "\n";
        }
    }

    ss <<
       "\n\tTotal unique contaminants: "                                     << count_contam       <<
          "(" << contam_percent << "%): "                                    <<
       "\n\t\tTranscriptome reference sequences labeled as a contaminant (FASTA):\n\t\t\t"<< out_best_contams_fa_prot<<
       "\n\t\tTranscriptome reference sequences labeled as a contaminant (TSV):\n\t\t\t"  << out_best_contams_tsv;


    // ********** Contaminant Calculations ************** //
    if (count_contam > 0) {
        ss << "\n\t\tFlagged contaminants (all % based on total contaminants):";
        for (auto &pair : contam_map) {
            percent = ((double)pair.second / count_contam) * 100;
            ss
                << "\n\t\t\t" << pair.first << ": " << pair.second << "(" << percent <<"%)";
        }
        ss << "\n\t\tTop 10 contaminants by species:";
        ct = 1;
        for (count_pair pair : contam_species_vect) {
            if (ct > 10) break;
            percent = ((double)pair.second / count_contam) * 100;
            ss
                << "\n\t\t\t" << ct << ")" << pair.first << ": "
                << pair.second << "(" << percent <<"%)";
            graph_contam_file << pair.first << '\t' << std::to_string(pair.second) << std::endl;
            ct++;
        }
    }

    ss << "\n\tTop 10 alignments by species:";
    ct = 1;
    for (count_pair pair : species_vect) {
        if (ct > 10) break;
        percent = ((double)pair.second / count_filtered) * 100;
        ss
            << "\n\t\t\t" << ct << ")" << pair.first << ": "
            << pair.second << "(" << percent <<"%)";
        graph_species_file << pair.first << '\t' << std::to_string(pair.second) << std::endl;
        ct++;
    }

    // ********* Graphing Handle ********** //
    graphingStruct.software_flag = GRAPH_SOFTWARE_FLAG;
    graph_contam_file.close();
    graph_species_file.close();
    graph_sum_file.close();
    if (count_contam > 0) {
        graphingStruct.fig_out_path   = graph_contam_png_path;
        graphingStruct.graph_title    = database + GRAPH_CONTAM_TITLE;
        graphingStruct.text_file_path = graph_contam_txt_path;
        graphingStruct.graph_type     = GRAPH_BAR_FLAG;
        _pGraphingManager->graph(graphingStruct);
    }
    graphingStruct.fig_out_path   = graph_species_png_path;
    graphingStruct.graph_title    = database + GRAPH_SPECIES_TITLE;
    graphingStruct.text_file_path = graph_species_txt_path;
    graphingStruct.graph_type     = GRAPH_BAR_FLAG;
    _pGraphingManager->graph(graphingStruct);

    graphingStruct.fig_out_path   = graph_sum_png_path;
    graphingStruct.graph_title    = database + GRAPH_DATABASE_SUM_TITLE;
    graphingStruct.text_file_path = graph_sum_txt_path;
    graphingStruct.graph_type     = GRAPH_SUM_FLAG;
    _pGraphingManager->graph(graphingStruct);

    // check if final - different graph
    // ************************************ ///

    return std::pair<std::string,std::string>(out_best_hits_fa_prot,out_no_hits_fa_prot);
}

/**
 *
 * @param diamond_maps
 * @return - pair of best_hit.fasta, no_hit.fasta
 */
std::pair<std::string,std::string> SimilaritySearch::process_best_diamond_hit(std::list<std::map<std::string,QuerySequence>> &diamond_maps) {
    print_debug("Compiling similarity results results to find best overall hits...");

    std::pair<std::string,std::string>  out_pair;
    std::stringstream                   out_stream;
    std::string                         out_msg;
    std::map<std::string,QuerySequence> compiled_hit_map;

    boostFS::remove_all(_results_path.c_str());
    boostFS::create_directories(_results_path.c_str());
    for (std::map<std::string,QuerySequence> &database_map : diamond_maps) {
        for (auto &pair : database_map) {
            pair.second.setIs_better_hit(false);
            pair.second.set_is_database_hit(true);
            std::map<std::string,QuerySequence>::iterator it = compiled_hit_map.find(pair.first);
            if (it != compiled_hit_map.end()) {
                if (pair.second > it->second) it->second = pair.second;
            } else {
                compiled_hit_map.emplace(pair);
            }
        }
    }
    out_stream<<std::fixed<<std::setprecision(2);
    out_stream << ENTAP_STATS::SOFTWARE_BREAK
               << "Compiled Similarity Search - Diamond - Best Overall\n"
               << ENTAP_STATS::SOFTWARE_BREAK;
    out_pair = calculate_best_stats(compiled_hit_map,out_stream,_results_path,true);
    out_msg  = out_stream.str() + "\n";
    print_statistics(out_msg);
    diamond_maps.clear();
    print_debug("Success!");
    return out_pair;
}

std::unordered_map<std::string, std::string> SimilaritySearch::read_tax_map() {
    print_debug("Reading taxonomic database into memory...");

    std::unordered_map<std::string, std::string> restored_map;

    if (!file_exists(TAX_DB_PATH)) {
        throw ExceptionHandler("NCBI Taxonomic database not found at: " +
            TAX_DB_PATH,ENTAP_ERR::E_INIT_TAX_READ);
    }
    try {
        {
            std::ifstream ifs(TAX_DB_PATH);
            boost::archive::binary_iarchive ia(ifs);
            ia >> restored_map;
        }
    } catch (std::exception &exception) {
        throw ExceptionHandler(exception.what(), ENTAP_ERR::E_INIT_TAX_READ);
    }
    print_debug("Success!");
    return restored_map;
}

std::pair<bool,std::string> SimilaritySearch::is_contaminant(std::string lineage, std::unordered_map<std::string, std::string> &database,
                    std::vector<std::string> &contams) {
    // species and tax database both lowercase
    if (contams.empty()) return std::pair<bool,std::string>(false,"");
    std::transform(lineage.begin(), lineage.end(), lineage.begin(), ::tolower);
    for (auto const &contaminant:contams) {
        if (lineage.find(contaminant) != std::string::npos){
            return std::pair<bool,std::string>(true,contaminant);
        }
    }
    return std::pair<bool,std::string>(false,"");
}

std::string SimilaritySearch::get_species(std::string &title) {
    // TODO use regex(database specific)

    std::string species;
    boost::smatch match;

    boost::regex ncbi_exp(_NCBI_REGEX);
    boost::regex uniprot_exp(_UNIPROT_REGEX);

    species = "";
    if (boost::regex_search(title,match,uniprot_exp)) {
        species = std::string(match[1].first, match[1].second);
    } else {
        if (boost::regex_search(title, match, ncbi_exp))
            species = std::string(match[1].first, match[1].second);
    }
    // Double bracket fix
    if (species[0] == '[') species = species.substr(1);
    if (species[species.length()-1] == ']') species = species.substr(0,species.length()-1);
    return species;
}

bool SimilaritySearch::is_informative(std::string title) {
    std::transform(title.begin(),title.end(),title.begin(),::tolower);
    for (std::string item : INFORMATIVENESS) {
        std::transform(item.begin(),item.end(),item.begin(),::tolower);
        if (title.find(item) != std::string::npos) return false;
    }
    return true;
}

void SimilaritySearch::print_header(std::string file) {
    std::ofstream ofstream(file, std::ios::out | std::ios::app);
    for (const std::string *header : DEFAULT_HEADERS) {
        ofstream << *header << "\t";
    }
    ofstream << std::endl;
    ofstream.close();
}

std::string SimilaritySearch::get_lineage(std::string species,
                                          std::unordered_map<std::string, std::string>&database) {
    if (species.empty()) return "";
    std::transform(species.begin(), species.end(), species.begin(), ::tolower);
    std::string lineage = "";
    if (database.find(species) != database.end() && !database.at(species).empty()) {
        lineage = database[species];
    } else {
        std::string temp_species = species;
        while (true) {
            unsigned long index = temp_species.find_last_of(" ");
            if (index == std::string::npos)break;
            temp_species = temp_species.substr(0,index);
            if (database.find(temp_species) != database.end()) {
                lineage = database[temp_species];
                break;
            }
        }
    }
    if (lineage.find("||") != std::string::npos) {
        return lineage.substr(lineage.find("||")+2);
    } else return lineage;
}