//
// Created by harta on 5/10/17.
//

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <csv.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/regex.hpp>
#include <iomanip>
#include "SimilaritySearch.h"
#include "ExceptionHandler.h"
#include "EntapInit.h"
#include "EntapConsts.h"
#include "QuerySequence.h"
#include "EntapExecute.h"



namespace boostFS = boost::filesystem;

SimilaritySearch::SimilaritySearch(std::vector<std::string> &databases, std::string input,
                           int threads, std::string exe, std::string out,std::string entap_exe,
                           boost::program_options::variables_map &user_flags, GraphingManager *graphingManager) {
    _database_paths = databases;
    _input_path = input;
    _threads = threads;
    _diamond_exe = exe;
    _outpath = out;
    _entap_exe = entap_exe;
    _input_species = "";
    if (user_flags.count(ENTAP_CONFIG::INPUT_FLAG_SPECIES)) {
        _input_species = user_flags[ENTAP_CONFIG::INPUT_FLAG_SPECIES].as<std::string>();
        std::transform(_input_species.begin(), _input_species.end(), _input_species.begin(), ::tolower);
        _input_species =
                _input_species.substr(0,_input_species.find("_")) +
                " " + _input_species.substr(_input_species.find("_")+1);
    }
    _qcoverage = user_flags[ENTAP_CONFIG::INPUT_FLAG_QCOVERAGE].as<double>();
    _tcoverage = user_flags[ENTAP_CONFIG::INPUT_FLAG_TCOVERAGE].as<double>();
    _overwrite = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _e_val = user_flags["e"].as<double>();

    std::vector<std::string> contaminants;
    if (user_flags.count("contam")) {
        contaminants = user_flags["contam"].as<std::vector<std::string>>();
    }
    for (int ind = 0; ind < contaminants.size(); ind++) {
        if (contaminants[ind].empty()) continue;
        std::string &str = contaminants[ind];
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    }
    _contaminants = contaminants;
    _software_flag = 0;

    _sim_search_path = (boostFS::path(out) / boostFS::path(SIM_SEARCH_DIR)).string();
    _processed_path  = (boostFS::path(_sim_search_path) / boostFS::path(PROCESSED_DIR)).string();
    _results_path    = (boostFS::path(_sim_search_path) / boostFS::path(RESULTS_DIR)).string();
    _figure_path     = (boostFS::path(_sim_search_path) / boostFS::path(FIGURE_DIR)).string();
    _graphingManager = graphingManager;
}

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

std::pair<std::string,std::string> SimilaritySearch::parse_files(std::string new_input,
                                   std::map<std::string, QuerySequence>& MAP) {
    _input_path = new_input;
    try {
        switch (_software_flag) {
            case 0:
                return diamond_parse(_contaminants,MAP);
            default:
                return diamond_parse(_contaminants,MAP);
        }
    } catch (ExceptionHandler &e) {throw e;}
}

std::vector<std::string> SimilaritySearch::diamond() {
    entapInit::print_msg("Beginning to execute DIAMOND...");
    std::vector<std::string> out_paths;
    if (!entapInit::file_exists(_input_path)) {
        throw ExceptionHandler("Transcriptome file not found",ENTAP_ERR::E_INIT_INDX_DATA_NOT_FOUND);
    }
    boostFS::path transc_name(_input_path); transc_name=transc_name.stem();
    if (transc_name.has_stem()) transc_name = transc_name.stem(); //.fasta.faa
    if (_overwrite) {
        boostFS::remove_all(_sim_search_path);
    }
    boostFS::create_directories(_sim_search_path);
    // database verification already ran, don't need to verify each path
    try {
        // assume all paths should be .dmnd
        for (std::string data_path : _database_paths) {
            entapInit::print_msg("Searching against database located at: " + data_path + "...");
            boostFS::path database_name(data_path); database_name=database_name.stem();
            std::string filename = _blast_type + "_" + transc_name.string() + "_" +
                                   database_name.string();
            std::string out_path = (boostFS::path(_sim_search_path) / filename).string() + ".out";
            std::string std_out  = (boostFS::path(_sim_search_path) / filename).string() + "_std";
            _file_to_database[out_path] = database_name.string();
            if (entapInit::file_exists(out_path)) {
                entapInit::print_msg("File found at " + out_path + " skipping execution against this database");
                out_paths.push_back(out_path);
                continue;
            }
            diamond_blast(_input_path, out_path, std_out,data_path, _threads, _blast_type);
            entapInit::print_msg("Success! Results written to " + out_path);
            out_paths.push_back(out_path);
        }
    } catch (const ExceptionHandler &e) {throw e;}
    _sim_search_paths = out_paths;
    return out_paths;
}

void SimilaritySearch::diamond_blast(std::string input_file, std::string output_file, std::string std_out,
                   std::string &database,int &threads, std::string &blast) {
    std::string diamond_run = _diamond_exe + " " + blast +" -d " + database + " --query-cover " + std::to_string(_qcoverage) +
                              " --subject-cover " + std::to_string(_tcoverage) +
                              " --more-sensitive" + " --top 3" + " -q " + input_file + " -o " + output_file + " -p " + std::to_string(threads) +" -f " +
                              "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle";
    entapInit::print_msg("\nExecuting Diamond:\n" + diamond_run);
    if (entapInit::execute_cmd(diamond_run, std_out) != 0) {
        throw ExceptionHandler("Error in DIAMOND run with database located at: " +
                               database, ENTAP_ERR::E_INIT_TAX_INDEX);
    }
}

std::vector<std::string> SimilaritySearch::verify_diamond_files(std::string &outpath, std::string name) {
    entapInit::print_msg("Override unselected, checking for diamond files"
                                 " of selected databases...");
    std::vector<std::string> out_list;
    for (std::string data_path : _database_paths) {
        // assume all paths should be .dmnd
        boostFS::path file_name(data_path);
        file_name = file_name.stem();
        std::string temp_out = outpath + _blast_type + "_" + name + "_" +
                               file_name.string() + ".out";
        if (!entapInit::file_exists(temp_out)){
            entapInit::print_msg("File at: " + temp_out + " not found, running diamond");
            out_list.clear();return out_list;
        }
        out_list.push_back(temp_out);
    }
    entapInit::print_msg("All diamond files found, skipping this stage of enTAP");
    _sim_search_paths = out_list;
    return out_list;
}

// input: 3 database string array of selected databases
std::pair<std::string,std::string> SimilaritySearch::diamond_parse(std::vector<std::string>& contams,
                                                                   std::map<std::string, QuerySequence> &SEQUENCES) {
    entapInit::print_msg("Beginning to filter individual diamond_files...");

    std::unordered_map<std::string, std::string>    taxonomic_database;
    std::list<std::map<std::string,QuerySequence>>  database_maps;

    try {
        taxonomic_database = read_tax_map();
    } catch (ExceptionHandler &e) {throw e;}
    if (_sim_search_paths.empty()) {
        std::string diamond_out = _outpath + ENTAP_CONFIG::SIM_SEARCH_OUT_PATH; // Raw outpath
        boostFS::path transc_name(_input_path); transc_name=transc_name.stem();
        if (transc_name.has_stem()) transc_name = transc_name.stem(); //.fasta.faa
        _sim_search_paths = verify_diamond_files(diamond_out,transc_name.string());
    }
    if (_sim_search_paths.empty()) throw ExceptionHandler("No diamond files found", ENTAP_ERR::E_RUN_SIM_SEARCH_FILTER);

    boostFS::remove_all(_processed_path);
    boostFS::create_directories(_processed_path);
    _input_lineage = get_lineage(_input_species,taxonomic_database);

    for (std::string &data : _sim_search_paths) {
        entapInit::print_msg("Diamond file located at " + data + " being filtered");
        io::CSVReader<ENTAP_EXECUTE::diamond_col_num, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(data);
        // todo have columns from input file, in_read_header for versatility
        std::string qseqid, sseqid, stitle, database_name;
        double evalue, pident, bitscore, coverage;
        int length, mismatch, gapopen, qstart, qend, sstart, send;
        unsigned long count_removed=0, count_TOTAL_hits=0, count_under_e=0;

        if (_file_to_database.find(data) != _file_to_database.end()) {
            database_name = _file_to_database[data];
        } else database_name = boostFS::path(data).filename().stem().string();
        std::string out_base_path = (boostFS::path(_processed_path) / boostFS::path(database_name)).string();
        std::string out_unselected_tsv = out_base_path + SIM_SEARCH_DATABASE_UNSELECTED;

        std::ofstream file_unselected_tsv(out_unselected_tsv,std::ios::out | std::ios::app);
        std::map<std::string, QuerySequence> database_map;
        print_header(out_unselected_tsv);

        while (in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                           qstart, qend, sstart, send, evalue, bitscore, coverage,stitle)) {
            count_TOTAL_hits++;

            QuerySequence new_query = QuerySequence();
            new_query.set_sim_search_results(data, qseqid, sseqid, pident, length, mismatch, gapopen,
                                             qstart, qend, sstart, send, evalue, bitscore, coverage, stitle);
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
            new_query.setFrame(SEQUENCES[qseqid].getFrame());  // May want to handle differently, SLOW
            new_query.set_tax_score(_input_lineage);
            if (evalue > _e_val) {
                count_under_e++; count_removed++;
                file_unselected_tsv << new_query <<std::endl;
                continue;
            }
            std::map<std::string, QuerySequence>::iterator it = database_map.find(qseqid);
            if (it != database_map.end()) {
                if (new_query > it->second) {
                    file_unselected_tsv << it->second << std::endl;
                    it->second = new_query;
                    count_removed++;
                } else {
                    count_removed++;
                    file_unselected_tsv << new_query<< std::endl;
                }
            } else database_map.emplace(qseqid, new_query);
        }
        entapInit::print_msg("File parsed, calculating statistics and writing output...");
        // final database stats
        file_unselected_tsv.close();
        std::stringstream out_stream;out_stream<<std::fixed<<std::setprecision(2);
        out_stream << ENTAP_STATS::SOFTWARE_BREAK
                   << "Similarity Search - Diamond - "<<data<<"\n"
                   << ENTAP_STATS::SOFTWARE_BREAK <<
                   "Statistics of file located at: "              << data               <<
                   "\n\tTotal hits: "                             << count_TOTAL_hits   <<
                   "\n\tUnselected results: "                     << count_removed      <<
                   "\n\t\tWritten to: "                           << out_unselected_tsv;
        calculate_best_stats(SEQUENCES,database_map,out_stream,out_base_path,false);
        std::string out_msg = out_stream.str() + "\n";
        entapExecute::print_statistics(out_msg,_outpath);
        database_maps.push_back(database_map);
        entapInit::print_msg("Success!");
    }
    return process_best_diamond_hit(database_maps,SEQUENCES);
}

std::pair<std::string,std::string> SimilaritySearch::calculate_best_stats (std::map<std::string, QuerySequence>&SEQUENCES,
                                             std::map<std::string, QuerySequence>&best_hits,
                                             std::stringstream &ss, std::string &base_path, bool is_final) {

    GraphingStruct              graphingStruct;
    std::string                 database;
    std::string                 figure_base;
    unsigned long               count_no_hit=0;
    unsigned long               count_contam=0;
    unsigned long               count_filtered=0;
    unsigned long               count_informative=0;
    unsigned long               count_uninformative=0;
    double                      contam_percent;
    std::map<std::string, int>  contam_map;
    std::map<std::string, int>  species_map;
    std::map<std::string, int>  contam_species_map;

    typedef std::pair<std::string,int> count_pair;
    struct compair {
        bool operator ()(count_pair const& one, count_pair const& two) const {
            return one.second > two.second;
        }
    };

    database = boostFS::path(base_path).filename().string();
    figure_base = (boostFS::path(_figure_path) / database).string();

    std::string out_best_contams_tsv             = base_path + SIM_SEARCH_DATABASE_CONTAM_TSV;
    std::string out_best_contams_fa_nucl         = base_path + SIM_SEARCH_DATABASE_CONTAM_FA_NUCL;
    std::string out_best_contams_fa_prot         = base_path + SIM_SEARCH_DATABASE_CONTAM_FA_PROT;

    std::string out_best_hits_tsv                = base_path + SIM_SEARCH_DATABASE_BEST_TSV;
    std::string out_best_hits_no_contam_tsv      = base_path + SIM_SEARCH_DATABASE_BEST_TSV_NO_CONTAM;
    std::string out_best_hits_fa_nucl            = base_path + SIM_SEARCH_DATABASE_BEST_FA_NUCL;
    std::string out_best_hits_fa_prot            = base_path + SIM_SEARCH_DATABASE_BEST_FA_PROT;
    std::string out_best_hits_fa_nucl_no_contam  = base_path + SIM_SEARCH_DATABASE_BEST_FA_NUCL_NO_CONTAM;
    std::string out_best_hits_fa_prot_no_contam  = base_path + SIM_SEARCH_DATABASE_BEST_FA_PROT_NO_CONTAM;

    std::string out_no_hits_fa_nucl              = base_path + SIM_SEARCH_DATABASE_NO_HITS_NUCL;
    std::string out_no_hits_fa_prot              = base_path + SIM_SEARCH_DATABASE_NO_HITS_PROT;

    std::string graph_species_txt_path           = figure_base + GRAPH_SPECIES_BAR_TXT;
    std::string graph_species_png_path           = figure_base + GRAPH_SPECIES_BAR_PNG;
    std::string graph_contam_txt_path            = figure_base + GRAPH_CONTAM_BAR_TXT;
    std::string graph_contam_png_path            = figure_base + GRAPH_CONTAM_BAR_PNG;

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

    graph_species_file << "Species\tCount"     << std::endl;
    graph_contam_file  << "Contaminant\tCount" << std::endl;

    for (auto &pair : SEQUENCES) {
        std::map<std::string, QuerySequence>::iterator it = best_hits.find(pair.first);
        if (it == best_hits.end()) {
            if (pair.second.isIs_protein()) {
                count_no_hit++;
                file_no_hits_nucl << pair.second.get_sequence_n() <<std::endl;
                file_no_hits_prot << pair.second.get_sequence_p() <<std::endl;
            }
        } else {
            file_best_hits_fa_nucl << pair.second.get_sequence_n()<<std::endl;
            file_best_hits_fa_prot << pair.second.get_sequence_p()<<std::endl;
            file_best_hits_tsv << it->second << std::endl;

            count_filtered++;
            std::string species;
            if (!it->second.get_species().empty()) {
                species = it->second.get_species();
            }
            if (it->second.isContaminant()) {
                count_contam++;
                file_best_contam_fa_nucl << pair.second.get_sequence_n()<<std::endl;
                file_best_contam_fa_prot << pair.second.get_sequence_p()<<std::endl;
                file_best_contam_tsv << it->second << std::endl;
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
                file_best_hits_tsv_no_contam << it->second << std::endl;
            }
            if (species_map.count(species)) {
                species_map[species]++;
            } else species_map[species] = 1;

            if (it->second.is_informative()) {
                count_informative++;
            } else count_uninformative++;

            if (is_final) {
                // TODO fix combining objects, move to class
                it->second.setSeq_length(pair.second.getSeq_length());
                it->second.set_sequence_n(pair.second.get_sequence_n());
                it->second.set_sequence_p(pair.second.get_sequence_p());
                it->second.setIs_protein(pair.second.isIs_protein());
                it->second.set_is_expression_kept(pair.second.is_is_expression_kept());
                pair.second = it->second;
            }
        }
    }
    file_best_hits_tsv.close();file_best_hits_tsv_no_contam.close();file_best_hits_fa_nucl.close();
    file_best_hits_fa_prot.close();file_best_hits_fa_nucl_no_contam.close();file_best_hits_fa_prot_no_contam.close();
    file_best_contam_tsv.close();file_best_contam_fa_prot.close();file_best_contam_fa_nucl.close();
    file_no_hits_nucl.close();file_no_hits_prot.close();

    std::vector<count_pair> contam_species_vect(contam_species_map.begin(),contam_species_map.end());
    std::vector<count_pair> species_vect(species_map.begin(),species_map.end());
    std::sort(contam_species_vect.begin(),contam_species_vect.end(),compair());
    std::sort(species_vect.begin(),species_vect.end(),compair());
    contam_percent = ((double)count_contam / count_filtered) * 100;

    ss <<
       "\n\tUnique hits: "                            << count_filtered     <<
       "\n\t\tBest fasta hits written to: "           << out_best_hits_fa_prot   <<
       "\n\t\tBest tsv hits written to: "             << out_best_hits_tsv  <<
       "\n\tSequences that did not hit: "             << count_no_hit       <<
       "\n\t\tWritten to: "                           << out_no_hits_fa_prot     <<
       "\n\tInformative hits: "                       << count_informative  <<
       "\n\tUninformative hits: "                     << count_uninformative<<
       "\n\tContaminants: "                           << count_contam       <<
          "(" << contam_percent << "): "                                    <<
       "\n\t\tFasta contaminants written to: "        << out_best_contams_fa_prot<<
       "\n\t\tTsv contaminants written to: "          << out_best_contams_tsv;

    if (count_contam > 0) {
        ss << "\n\t\tFlagged contaminants (all % based on total contaminants):";
        for (auto &pair : contam_map) {
            double percent = ((double)pair.second / count_contam) * 100;
            ss
                << "\n\t\t\t" << pair.first << ": " << pair.second << "(" << percent <<"%)";
        }
        ss << "\n\t\tTop 10 contaminants by species:";
        int ct = 1;
        for (count_pair pair : contam_species_vect) {
            if (ct > 10) break;
            double percent = ((double)pair.second / count_contam) * 100;
            ss
                << "\n\t\t\t" << ct << ")" << pair.first << ": "
                << pair.second << "(" << percent <<"%)";
            graph_contam_file << pair.first << '\t' << std::to_string(pair.second) << std::endl;
            ct++;
        }
    }

    ss << "\n\tTop 10 species:";
    int ct = 1;
    for (count_pair pair : species_vect) {
        if (ct > 10) break;
        double percent = ((double)pair.second / count_filtered) * 100;
        ss
            << "\n\t\t\t" << ct << ")" << pair.first << ": "
            << pair.second << "(" << percent <<"%)";
        graph_species_file << pair.first << '\t' << std::to_string(pair.second) << std::endl;
        ct++;
    }

    // Graphing Handle
    graphingStruct.software_flag = GRAPH_SOFTWARE_FLAG;
    graph_contam_file.close();
    graph_species_file.close();
    if (is_final) {

    } else {
        if (count_contam > 0) {
            graphingStruct.fig_out_path = graph_contam_png_path;
            graphingStruct.graph_title  = database + GRAPH_CONTAM_TITLE;
            graphingStruct.text_file_path = graph_contam_txt_path;
            graphingStruct.graph_type = GRAPH_BAR_FLAG;
            _graphingManager->graph(graphingStruct);
        }
        graphingStruct.fig_out_path = graph_species_png_path;
        graphingStruct.graph_title  = database + GRAPH_SPECIES_TITLE;
        graphingStruct.text_file_path = graph_species_txt_path;
        graphingStruct.graph_type = GRAPH_BAR_FLAG;
        _graphingManager->graph(graphingStruct);
    }





    return std::pair<std::string,std::string>(out_best_hits_fa_prot,out_no_hits_fa_prot);
}

/**
 *
 * @param diamond_maps
 * @return - pair of best_hit.fasta, no_hit.fasta
 */
std::pair<std::string,std::string> SimilaritySearch::process_best_diamond_hit(std::list<std::map<std::string,QuerySequence>> &diamond_maps,
                                      std::map<std::string, QuerySequence>&SEQUENCES) {
    entapInit::print_msg("Compiling similarity results results to find best overall hits...");
    boostFS::remove_all(_results_path.c_str());
    boostFS::create_directories(_results_path.c_str());
    std::map<std::string,QuerySequence> compiled_hit_map;
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
    std::stringstream out_stream;out_stream<<std::fixed<<std::setprecision(2);
    out_stream <<
            "------Compiled Results (Best hit selection across all databases ------";
    std::pair<std::string,std::string> out_pair =
            calculate_best_stats(SEQUENCES,compiled_hit_map,out_stream,_results_path,true);
    std::string out_msg = out_stream.str() + "\n";
    entapExecute::print_statistics(out_msg,_outpath);
    diamond_maps.clear();
    entapInit::print_msg("Success!");
    return out_pair;
}

std::unordered_map<std::string, std::string> SimilaritySearch::read_tax_map() {
    entapInit::print_msg("Reading taxonomic database into memory...");
    std::unordered_map<std::string, std::string> restored_map;
    std::string tax_path = _entap_exe + ENTAP_CONFIG::TAX_BIN_PATH;
    if (!entapInit::file_exists(tax_path)) {
        throw ExceptionHandler("NCBI Taxonomic database not found at: " +
            tax_path,ENTAP_ERR::E_INIT_TAX_READ);
    }
    try {
        {
            std::ifstream ifs(tax_path);
            boost::archive::binary_iarchive ia(ifs);
            ia >> restored_map;
        }
    } catch (std::exception &exception) {
        throw ExceptionHandler(exception.what(), ENTAP_ERR::E_INIT_TAX_READ);
    }
    entapInit::print_msg("Success!");
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
    boost::regex ncbi_exp(_NCBI_REGEX);
    boost::regex uniprot_exp(_UNIPROT_REGEX);
    boost::smatch match;
    std::string species = "";
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
    ofstream <<
         "Query Seq\t"
                 "Subject Seq\t"
                 "Percent Identical\t"
                 "Alignment Length\t"
                 "Mismatches\t"
                 "Gap Openings\t"
                 "Query Start\t"
                 "Query End\t"
                 "Subject Start\t"
                 "Subject Eng\t"
                 "E Value\t"
                 "Coverage\t"
                 "Informativeness\t"
                 "Species\t"
                 "Origin Database\t"
                 "Frame\t"

         <<std::endl;
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