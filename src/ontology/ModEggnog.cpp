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


#include <boost/archive/binary_iarchive.hpp>
#include <iomanip>
#include <fstream>
#include <csv.h>
#include "ModEggnog.h"
#include "../ExceptionHandler.h"

std::pair<bool, std::string> ModEggnog::verify_files() {
    std::string                        annotation_base_flag;
    std::string                        annotation_no_flag;
    bool                               verified;

    annotation_base_flag = PATHS(_egg_out_dir, "annotation_results");
    annotation_no_flag   = PATHS(_egg_out_dir, "annotation_results_no_hits");
    _out_hits            = annotation_base_flag  +".emapper.annotations";
    _out_no_hits         = annotation_no_flag +".emapper.annotations";

    verified = false;
    print_debug("Overwrite was unselected, verifying output files...");
    if (file_exists(_out_hits)) {
        print_debug("File located at: " + _out_hits + " found");
        verified = true;
    } else print_debug("File located at: " + _out_hits + " NOT found");
    if (file_exists(_out_no_hits)) {
        print_debug("File located at: " + _out_no_hits + " found");
        verified = true;
    } else print_debug("File located at: " + _out_no_hits + " NOT found");
    if (verified) {
        print_debug("One or more ontology files were found, skipping ontology execution");
        return std::make_pair(true, "");
    } else {
        print_debug("No ontology files were found, continuing with execution");
        return std::make_pair(false,"");
    }
}


/**
 *
 * @param SEQUENCES
 * @param out - <egg results of sequences that hit databases, results of no hits>
 */
void ModEggnog::parse() {

    print_debug("Beginning to parse eggnog results...");

    typedef std::map<std::string,std::map<std::string, uint32>> GO_top_map_t;

    std::stringstream                        ss;
    std::string                              out_msg;
    std::string                              out_no_hits_nucl;
    std::string                              out_no_hits_prot;
    std::string                              out_hit_nucl;
    std::string                              out_hit_prot;
    std::string                              path;
    std::string                              fig_txt_bar_go_overall;
    std::string                              fig_png_bar_go_overall;
    std::string                              fig_txt_bar_ortho;
    std::string                              fig_png_bar_ortho;
    std::string                              tax_scope_readable;
    std::string                              fig_txt_go_bar;
    std::string                              fig_png_go_bar;
    std::map<std::string, struct_go_term>    GO_DATABASE;
    std::map<std::string, int>               eggnog_map;
    uint32                                   count_total_go_hits=0;
    uint32                                   count_total_go_terms=0;
    uint32                                   count_go_bio=0;
    uint32                                   count_go_cell=0;
    uint32                                   count_go_mole=0;
    uint32                                   count_no_go=0;
    uint32                                   count_no_kegg=0;
    uint32                                   count_TOTAL_hits=0;         // All ortho matches
    uint32                                   count_total_kegg_terms=0;
    uint32                                   count_total_kegg_hits=0;
    uint32                                   count_no_hits=0;            // Unannotated OGs
    uint32                                   count_tax_scope=0;
    uint32                                   ct = 0;
    fp32                                     percent;
    DatabaseHelper                           EGGNOG_DATABASE;
    std::map<std::string, uint32>            tax_scope_ct_map;
    GO_top_map_t                             go_combined_map;     // Just for convenience
    go_struct                                go_parsed;
    GraphingStruct                           graphingStruct;

    ss<<std::fixed<<std::setprecision(2);
    boostFS::remove_all(_proc_dir);
    boostFS::create_directories(_proc_dir);
    boostFS::create_directories(_figure_dir);
    try {
        GO_DATABASE = read_go_map();
    } catch (ExceptionHandler const &e) {throw e;}

    if (!EGGNOG_DATABASE.open(_eggnog_db_path))
        throw ExceptionHandler("Unable to open GO database",ENTAP_ERR::E_PARSE_EGGNOG);

    for (int i=0; i<2;i++) {
        i == 0 ? path=_out_hits : path=_out_no_hits;
        print_debug("Eggnog file located at " + path + " being filtered");
        if (!file_exists(path)) {
            print_debug("File not found, skipping...");continue;
        }
        path = eggnog_format(path);
        std::string qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                best_og, cog_cat, eggnog_annot;
        io::CSVReader<EGGNOG_COL_NUM, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(path);
        while (in.read_row(qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                           best_og, cog_cat, eggnog_annot)) {
            query_map_struct::iterator it = (*pQUERY_DATA->get_sequences_ptr()).find(qseqid);
            if (it != (*pQUERY_DATA->get_sequences_ptr()).end()) {
                count_TOTAL_hits++;
                it->second.set_eggnog_results(seed_ortho,seed_e,seed_score,predicted_gene,go_terms,
                                              kegg,tax_scope,ogs, EGGNOG_DATABASE);
                go_parsed = parse_go_list(go_terms,GO_DATABASE,',');
                it->second.set_go_parsed(go_parsed, ENTAP_EXECUTE::EGGNOG_INT_FLAG);
                it->second.set_is_family_assigned(true);
                eggnog_map[qseqid] = 1;
                if (!go_parsed.empty()) {
                    count_total_go_hits++;
                    it->second.set_is_one_go(true);
                    for (auto &pair : go_parsed) {
                        for (std::string &term : pair.second) {
                            count_total_go_terms++;
                            if (pair.first == GO_MOLECULAR_FLAG) {
                                count_go_mole++;
                            } else if (pair.first == GO_CELLULAR_FLAG) {
                                count_go_cell++;
                            } else if (pair.first == GO_BIOLOGICAL_FLAG) {
                                count_go_bio++;
                            }
                            if (go_combined_map[pair.first].count(term)) {
                                go_combined_map[pair.first][term]++;
                            } else go_combined_map[pair.first][term] = 1;
                            if (go_combined_map[GO_OVERALL_FLAG].count(term)) {
                                go_combined_map[GO_OVERALL_FLAG][term]++;
                            } else go_combined_map[GO_OVERALL_FLAG][term]=1;
                        }
                    }
                } else {
                    count_no_go++;
                }

                // Compile KEGG information
                if (!kegg.empty()) {
                    count_total_kegg_hits++;
                    ct = (uint32) std::count(kegg.begin(), kegg.end(), ',');
                    count_total_kegg_terms += ct + 1;
                    it->second.set_is_one_kegg(true);
                } else {
                    count_no_kegg++;
                }

                // Compile Taxonomic Orthogroup stats
                tax_scope_readable = it->second.get_tax_scope();
                if (!tax_scope_readable.empty()) {
                    count_tax_scope++;
                    if (tax_scope_ct_map.count(tax_scope_readable)) {
                        tax_scope_ct_map[tax_scope_readable]++;
                    } else tax_scope_ct_map[tax_scope_readable] = 1;

                }
            }
        }
        boostFS::remove(path);
    }

    EGGNOG_DATABASE.close();
    out_no_hits_nucl = (boostFS::path(_proc_dir) / boostFS::path(OUT_UNANNOTATED_NUCL)).string();
    out_no_hits_prot = (boostFS::path(_proc_dir) / boostFS::path(OUT_UNANNOTATED_PROT)).string();
    out_hit_nucl     = (boostFS::path(_proc_dir) / boostFS::path(OUT_ANNOTATED_NUCL)).string();
    out_hit_prot     = (boostFS::path(_proc_dir) / boostFS::path(OUT_ANNOTATED_PROT)).string();
    std::ofstream file_no_hits_nucl(out_no_hits_nucl, std::ios::out | std::ios::app);
    std::ofstream file_no_hits_prot(out_no_hits_prot, std::ios::out | std::ios::app);
    std::ofstream file_hits_nucl(out_hit_nucl, std::ios::out | std::ios::app);
    std::ofstream file_hits_prot(out_hit_prot, std::ios::out | std::ios::app);

    print_debug("Success! Computing overall statistics...");
    for (auto &pair : *pQUERY_DATA->get_sequences_ptr()) {
        if (eggnog_map.find(pair.first) == eggnog_map.end()) {
            // Unannotated sequence
            if (!pair.second.get_sequence_n().empty()) file_no_hits_nucl<<pair.second.get_sequence_n()<<std::endl;
            file_no_hits_prot << pair.second.get_sequence_p() << std::endl;
            count_no_hits++;
        } else {
            // Annotated sequence
            if (!pair.second.get_sequence_n().empty()) file_hits_nucl<<pair.second.get_sequence_n()<<std::endl;
            file_hits_prot << pair.second.get_sequence_p() << std::endl;
        }
    }

    file_hits_nucl.close();
    file_hits_prot.close();
    file_no_hits_nucl.close();
    file_no_hits_prot.close();

    ss << ENTAP_STATS::SOFTWARE_BREAK                             <<
          "Gene Family - Gene Ontology and Pathway - Eggnog\n"    <<
          ENTAP_STATS::SOFTWARE_BREAK                             <<
          "Statistics for overall Eggnog results: "               <<
          "\nTotal unique sequences with family assignment: "     << count_TOTAL_hits <<
          "\nTotal unique sequences without family assignment: "  << count_no_hits;

    // -------- Top Ten Taxonomic Scopes ------- //
    if (!tax_scope_ct_map.empty()) {
        std::string fig_txt_tax_bar = (boostFS::path(_figure_dir) / GRAPH_EGG_TAX_BAR_TXT).string();
        std::string fig_png_tax_bar = (boostFS::path(_figure_dir) / GRAPH_EGG_TAX_BAR_PNG).string();
        std::ofstream file_tax_bar(fig_txt_tax_bar, std::ios::out | std::ios::app);
        file_tax_bar << "Taxonomic Scope\tCount" << std::endl;

        ss << "\nTop 10 Taxonomic Scopes Assigned:";
        ct = 1;
        std::vector<count_pair> tax_scope_vect(tax_scope_ct_map.begin(),tax_scope_ct_map.end());
        std::sort(tax_scope_vect.begin(),tax_scope_vect.end(),compair());
        for (count_pair &pair : tax_scope_vect) {
            if (ct > 10) break;
            percent = ((fp32)pair.second / count_tax_scope) * 100;
            ss <<
               "\n\t" << ct << ")" << pair.first << ": " << pair.second <<
               "(" << percent << "%)";
            file_tax_bar << pair.first << '\t' << std::to_string(pair.second) << std::endl;
            ct++;
        }
        file_tax_bar.close();
        graphingStruct.fig_out_path = fig_png_tax_bar;
        graphingStruct.text_file_path = fig_txt_tax_bar;
        graphingStruct.graph_title = GRAPH_EGG_TAX_BAR_TITLE;
        graphingStruct.software_flag = GRAPH_ONTOLOGY_FLAG;
        graphingStruct.graph_type = GRAPH_TOP_BAR_FLAG;
        pGraphingManager->graph(graphingStruct);
    }
    // --------------------------------------- //

    ss<<
      "\nTotal unique sequences with at least one GO term: " << count_total_go_hits <<
      "\nTotal unique sequences without GO terms: " << count_no_go <<
      "\nTotal GO terms assigned: " << count_total_go_terms;

    if (count_total_go_hits > 0) {
        for (uint16 lvl : _go_levels) {
            for (auto &pair : go_combined_map) {
                if (pair.first.empty()) continue;
                // Count maps (biological/molecular/cellular/overall)
                fig_txt_go_bar = (boostFS::path(_figure_dir) / pair.first).string() + std::to_string(lvl)+GRAPH_GO_END_TXT;
                fig_png_go_bar = (boostFS::path(_figure_dir) / pair.first).string() + std::to_string(lvl)+GRAPH_GO_END_PNG;
                std::ofstream file_go_bar(fig_txt_go_bar, std::ios::out | std::ios::app);
                std::vector<count_pair> go_vect(pair.second.begin(),pair.second.end());
                std::sort(go_vect.begin(),go_vect.end(),compair());
                file_go_bar << "Gene Ontology Term\tCount" << std::endl;

                // get total count for each level...change, didn't feel like making another
                uint32 lvl_ct = 0;   // Use for percentages, total terms for each lvl
                ct = 0;              // Use for unique count
                for (count_pair &pair2 : go_vect) {
                    if (pair2.first.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                        ct++;
                        lvl_ct += pair2.second;
                    }
                }
                ss << "\nTotal "        << pair.first <<" terms (lvl="          << lvl << "): " << lvl_ct;
                ss << "\nTotal unique " << pair.first <<" terms (lvl="          << lvl << "): " << ct;
                ss << "\nTop 10 "       << pair.first <<" terms assigned (lvl=" << lvl << "): ";

                ct = 1;
                for (count_pair &pair2 : go_vect) {
                    if (ct > 10) break;
                    if (pair2.first.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                        percent = ((fp32)pair2.second / lvl_ct) * 100;
                        ss <<
                           "\n\t" << ct << ")" << pair2.first << ": " << pair2.second <<
                           "(" << percent << "%)";
                        file_go_bar << pair2.first << '\t' << std::to_string(pair2.second) << std::endl;
                        ct++;
                    }
                }
                file_go_bar.close();
                graphingStruct.fig_out_path   = fig_png_go_bar;
                graphingStruct.text_file_path = fig_txt_go_bar;
                if (pair.first == GO_BIOLOGICAL_FLAG) graphingStruct.graph_title = GRAPH_GO_BAR_BIO_TITLE + "_Level:_"+std::to_string(lvl);
                if (pair.first == GO_CELLULAR_FLAG) graphingStruct.graph_title = GRAPH_GO_BAR_CELL_TITLE+ "_Level:_"+std::to_string(lvl);
                if (pair.first == GO_MOLECULAR_FLAG) graphingStruct.graph_title = GRAPH_GO_BAR_MOLE_TITLE+ "_Level:_"+std::to_string(lvl);
                if (pair.first == GO_OVERALL_FLAG) graphingStruct.graph_title = GRAPH_GO_BAR_ALL_TITLE+ "_Level:_"+std::to_string(lvl);
                // Other params can stay the same
                pGraphingManager->graph(graphingStruct);
            }
        }

    }
    ss<<
      "\nTotal unique sequences with at least one pathway (KEGG) assignment: " << count_total_kegg_hits<<
      "\nTotal unique sequences without pathways (KEGG): " << count_no_kegg<<
      "\nTotal pathways (KEGG) assigned: " << count_total_kegg_terms;
    out_msg = ss.str();
    print_statistics(out_msg);
    GO_DATABASE.clear();
    print_debug("Success!");
}

void ModEggnog::execute() {
    print_debug("Running eggnog...");

    std::string                        annotation_base_flag;
    std::string                        annotation_no_flag;
    std::string                        annotation_std;
    std::string                        eggnog_command;
    std::string                        hit_out;
    std::string                        no_hit_out;

    annotation_base_flag = PATHS(_egg_out_dir, "annotation_results");
    annotation_no_flag   = PATHS(_egg_out_dir, "annotation_results_no_hits");
    annotation_std       = PATHS(_egg_out_dir, "annotation_std");
    eggnog_command       = "python " + _exe_path + " ";

    std::unordered_map<std::string,std::string> eggnog_command_map = {
            {"-i",_inpath},
            {"--output",annotation_base_flag},
            {"--cpu",std::to_string(_threads)},
            {"-m", "diamond"}
    };
    if (!_blastp) eggnog_command_map["--translate"] = " ";
    if (file_exists(_inpath)) {
        for (auto &pair : eggnog_command_map)eggnog_command += pair.first + " " + pair.second + " ";
        print_debug("\nExecuting eggnog mapper against protein sequences that hit databases...\n"
                    + eggnog_command);
        if (execute_cmd(eggnog_command, annotation_std) !=0) {
            throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
        }
        print_debug("Success! Results written to: " + annotation_base_flag);
    } else {
        throw ExceptionHandler("No input file found at: " + _inpath,
                               ENTAP_ERR::E_RUN_EGGNOG);
    }
    if (file_exists(_in_no_hits)) {
        std::ifstream inFile(_in_no_hits);
        long line_num = std::count(std::istreambuf_iterator<char>(inFile),
                                   std::istreambuf_iterator<char>(), '\n');
        inFile.close();
        if (line_num >1) {
            eggnog_command_map["-i"] = _in_no_hits;
            eggnog_command_map["--output"] = annotation_no_flag;
            eggnog_command = "python " + _exe_path + " ";
            for (auto &pair : eggnog_command_map) eggnog_command += pair.first + " " + pair.second + " ";
            print_debug("\nExecuting eggnog mapper against protein sequences that did not hit databases...\n"
                        + eggnog_command);
            if (execute_cmd(eggnog_command, annotation_std) !=0) {
                throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
            }
        }
    }
    print_debug("Success!");
}

void ModEggnog::set_data(std::string & eggnog_databse, std::vector<std::string>&) {

    _eggnog_db_path = eggnog_databse;

    _egg_out_dir= PATHS(_ontology_dir, EGGNOG_DIRECTORY);
    _figure_dir = PATHS(_egg_out_dir, FIGURE_DIR);
    _proc_dir   = PATHS(_egg_out_dir, PROCESSED_OUT_DIR);

    boostFS::remove_all(_figure_dir);
    boostFS::remove_all(_proc_dir);

    boostFS::create_directories(_egg_out_dir);
    boostFS::create_directories(_figure_dir);
    boostFS::create_directories(_proc_dir);
}


// TODO remove
std::string ModEggnog::eggnog_format(std::string file) {

    std::string out_path;
    std::string line;

    out_path = file + "_alt";
    boostFS::remove(out_path);
    std::ifstream in_file(file);
    std::ofstream out_file(out_path);
    while (getline(in_file,line)) {
        if (line.at(0) == '#' || line.empty()) continue;
        out_file << line << std::endl;
    }
    in_file.close();
    out_file.close();
    return out_path;
}
