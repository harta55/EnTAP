/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
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

#include "ModEggnogDMND.h"
#include "../database/EggnogDatabase.h"
#include "../TerminalCommands.h"

ModEggnogDMND::ModEggnogDMND(std::string &ont_out, std::string &in_hits,
                             EntapDataPtrs &entap_data, std::string &exe, std::string sql_db_path)
        : AbstractOntology(in_hits, ont_out, entap_data, "EggNOG_DMND", exe) {
    FS_dprint("Spawn Object - ModEggnogDMND");

    _eggnog_db_path = sql_db_path;
    _software_flag = ENTAP_EXECUTE::EGGNOG_DMND_INT_FLAG;

}

std::pair<bool, std::string> ModEggnogDMND::verify_files() {
    bool file_exists = false;
    uint16 file_status = 0;

    FS_dprint("Overwrite was unselected, verifying output files...");
    _out_hits = get_output_dmnd_filepath(true);
    file_status = _pFileSystem->get_file_status(_out_hits);

    if (file_status != 0) {
        FS_dprint(_pFileSystem->print_file_status(file_status,_out_hits));
        FS_dprint("Errors in opening file, continuing with execution...");
    } else {
        file_exists = true;
    }
    return std::make_pair(file_exists, _out_hits);
}

void ModEggnogDMND::execute() {
    std::string                        std_out;
    std::string                        cmd;
    std::string                        blast;
    TerminalData                       terminalData;

    FS_dprint("Running EggNOG against Diamond database...");

    // Ensure both input path and EggNOG DMND database exist before continuing
    if (!_pFileSystem->file_exists(EGG_DMND_PATH)) {
        throw ExceptionHandler("EggNOG DIAMOND database not found at: " + EGG_DMND_PATH,
                               ERR_ENTAP_EGGNOG_FILES);
    }
    if (!_pFileSystem->file_exists(_in_hits)) {
        throw ExceptionHandler("Input transcriptome not found at: " + _in_hits, ERR_ENTAP_EGGNOG_FILES);
    }

    // Generate paths for DIAMOND run (out_hits set previously)
    std_out = get_output_dmnd_filepath(false) + "_" + FileSystem::EXT_STD;
    _blastp ? blast = "blastp" : blast = "blastx";

    //Run DIAMOND
    cmd =
            DIAMOND_EXE + " " +
            blast +
            " -d " + EGG_DMND_PATH +
            " --top 1"             +
            " --more-sensitive"    +
            " -q "                 + _in_hits   +
            " -o "                 + _out_hits  +
            " -p "                 + std::to_string(_threads) +
            " -f " + "6 qseqid sseqid pident length mismatch gapopen "
                    "qstart qend sstart send evalue bitscore qcovhsp stitle";

    terminalData.command        = cmd;
    terminalData.print_files    = true;
    terminalData.base_std_path  = std_out;

    if (TC_execute_cmd(terminalData) != 0) {
        // Error in run
        _pFileSystem->delete_file(_out_hits);
        FS_dprint("DIAMOND STD OUT:\n" + terminalData.out_stream.str());
        throw ExceptionHandler("Error in running DIAMOND against EggNOG database at: " +
                               EGG_DMND_PATH + "\nDIAMOND Error:\n" + terminalData.err_stream.str(), ERR_ENTAP_RUN_EGGNOG);
    }
}

void ModEggnogDMND::parse() {
    uint16         file_status=0;
    uint64         sequence_ct=0;   // dprintf sequence count
    std::stringstream stats_stream;

    FS_dprint("Parsing EggNOG DMND file located at: " + _out_hits);

    // Ensure file is valid
    file_status = _pFileSystem->get_file_status(_out_hits);
    if (file_status != 0) {
        throw ExceptionHandler(_pFileSystem->print_file_status(file_status,_out_hits),
                               ERR_ENTAP_PARSE_EGGNOG_DMND);
    }

    // File valid, continue
    FS_dprint("Beginning to parse EggNOG results...");
    ENTAP_STATS::ES_format_stat_stream(stats_stream, "Gene Family - Gene Ontology and Pathway - EggNOG");

#ifdef USE_FAST_CSV
    // ------------------ Read from DIAMOND output ---------------------- //
    std::string qseqid;
    std::string sseqid, stitle, database_name,pident, bitscore,
            length, mismatch, gapopen, qstart, qend, sstart, send;
    fp64 evalue;
    fp64 coverage;
    QuerySequence::EggnogResults eggnogResults;
    QuerySequence *querySequence;
    // ----------------------------------------------------------------- //
    // Begin using CSVReader lib to parse data
    try {
        io::CSVReader<DMND_COL_NUMBER, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(_out_hits);
        while (in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                           qstart, qend, sstart, send, evalue, bitscore, coverage,stitle)) {
            // Currently throwing away most DIAMOND results

            // Print progress to debug
            if (++sequence_ct % STATUS_UPDATE_HITS == 0) {
                FS_dprint("Alignments parsed: " + std::to_string(sequence_ct));
            }

            // Ensure we recognize the query sequence before continuing
            querySequence = _pQUERY_DATA->get_sequence(qseqid);
            if (querySequence == nullptr) {
                throw ExceptionHandler("Unable to find sequence " + qseqid + " in input transcriptome",
                                       ERR_ENTAP_PARSE_EGGNOG_DMND);
            }

            // Populate seed data from diamond
            eggnogResults = {};
            eggnogResults.seed_evalue = evalue;
            eggnogResults.seed_score  = bitscore;
            eggnogResults.seed_coverage = coverage;
            eggnogResults.seed_ortholog = sseqid;
            querySequence->add_alignment(GENE_ONTOLOGY, _software_flag, eggnogResults, EGG_DMND_PATH);

        } // End WHILE in.read_row

        if (sequence_ct > 0) {
            FS_dprint("Success!");
            calculate_stats(stats_stream);
        } else {
            // NO alignments against EggNOG !!!
            FS_dprint("WARNING: NO alignments against EggNOG!");
            stats_stream << "Warning: No alignments against EggNOG database" << std::endl;
        }

    } catch (const ExceptionHandler &e) {
        throw e;
    } catch (const std::exception &e) {
        throw ExceptionHandler(e.what(), ERR_ENTAP_PARSE_EGGNOG_DMND);
    }
#endif
}



void ModEggnogDMND::calculate_stats(std::stringstream &stream) {
    FS_dprint("Success! Calculating statistics and accessing database...");

    EggnogDatabase                     *eggnogDatabase;
    QuerySequence::EggnogResults       *eggnog_results;
    QuerySequence::EggnogDmndAlignment *best_hit;
    Compair<std::string>                                  tax_scope_counter;
    std::unordered_map<std::string,Compair<std::string>>  go_combined_map;
    GraphingData                        graphingStruct;

    uint64         ct_alignments=0;
    uint64         ct_no_alignment=0;
    uint64         ct_total_go_hits=0;      // Sequences that had atleast one go
    uint64         ct_total_kegg_hits=0;    // Sequences that had atleast one kegg
    uint64         ct_total_kegg_terms=0;
    uint32         ct = 0;
    fp32           percent;

    std::string    out_msg;


    // Output files
    std::string out_no_hits_nucl = PATHS(_proc_dir, OUT_UNANNOTATED_NUCL);
    std::string out_no_hits_prot = PATHS(_proc_dir, OUT_UNANNOTATED_PROT);
    std::string out_hit_nucl     = PATHS(_proc_dir, OUT_ANNOTATED_NUCL);
    std::string out_hit_prot     = PATHS(_proc_dir, OUT_ANNOTATED_PROT);
    std::ofstream file_no_hits_nucl(out_no_hits_nucl, std::ios::out | std::ios::app);
    std::ofstream file_no_hits_prot(out_no_hits_prot, std::ios::out | std::ios::app);
    std::ofstream file_hits_nucl(out_hit_nucl, std::ios::out | std::ios::app);
    std::ofstream file_hits_prot(out_hit_prot, std::ios::out | std::ios::app);


    // Generate EggNOG database
    eggnogDatabase = new EggnogDatabase(_pFileSystem, _pEntapDatabase);
    if (eggnogDatabase->open_sql(EGG_SQL_DB_PATH) != EggnogDatabase::ERR_EGG_OK) {
        throw ExceptionHandler("Unable to open EggNOG SQL Database", ERR_ENTAP_PARSE_EGGNOG_DMND);
    }


    // Parse through all query sequences
    for (auto &pair : *_pQUERY_DATA->get_sequences_ptr()) {
        // Check if each sequence is an eggnog alignment
        if (pair.second->hit_database(GENE_ONTOLOGY, _software_flag, EGG_DMND_PATH)) {
            // Yes, hit EggNOG database
            ct_alignments++;

            if (!pair.second->get_sequence_n().empty()) file_hits_nucl << pair.second->get_sequence_n() << std::endl;
            if (!pair.second->get_sequence_p().empty()) file_hits_prot << pair.second->get_sequence_p() << std::endl;

            best_hit = pair.second->get_best_hit_alignment<QuerySequence::EggnogDmndAlignment>
                    (GENE_ONTOLOGY, _software_flag, EGG_DMND_PATH);
            eggnog_results = best_hit->get_results();
            // get SQL data
            eggnogDatabase->get_eggnog_entry(*eggnog_results);

            //  Analyze Gene Ontology Stats
            if (!eggnog_results->parsed_go.empty()) {
                ct_total_go_hits++;
                for (auto &pair2 : eggnog_results->parsed_go) {
                    // pair - first: GO category, second; vector of terms
                    for (std::string &term : pair2.second) {
                        // Count the terms we've found for individual category
                        go_combined_map[pair2.first].add_value(term);
                        // Count the terms we've found overall (not category specific)
                        go_combined_map[GO_OVERALL_FLAG].add_value(term);
                    }
                }
            }

            // Analyze KEGG stats
            if (!eggnog_results->kegg.empty()) {
                ct_total_kegg_hits++;
                ct_total_kegg_terms += (uint32) STR_COUNT(eggnog_results->kegg, ',') + 1;
            }

            // Compile Taxonomic Orthogroup stats
            if (!eggnog_results->tax_scope_readable.empty()) {
                // Count number of unique taxonomic groups
                tax_scope_counter.add_value(eggnog_results->tax_scope_readable);
            }


        } else {
            // No, did not hit database
            ct_no_alignment++;
            if (!pair.second->get_sequence_n().empty()) file_no_hits_nucl << pair.second->get_sequence_n() << std::endl;
            if (!pair.second->get_sequence_p().empty()) file_no_hits_prot << pair.second->get_sequence_p() << std::endl;
        }
    } // END FOR LOOP

    // Close EggNOG database, files
    delete eggnogDatabase;
    file_hits_nucl.close();
    file_hits_prot.close();
    file_no_hits_nucl.close();
    file_no_hits_prot.close();

    // Begin to print stats / files
    stream <<
       "Statistics for overall Eggnog results: "               <<
       "\nTotal unique sequences with family assignment: "     << ct_alignments <<
       "\nTotal unique sequences without family assignment: "  << ct_no_alignment;

    // Make sure we have hits before doing anything

    //--------------------- Top Ten Taxonomic Scopes --------------//
    if (!tax_scope_counter.empty()) {
        // Setup graphing files
        std::string fig_txt_tax_bar = PATHS(_figure_dir, GRAPH_EGG_TAX_BAR_TXT);
        std::string fig_png_tax_bar = PATHS(_figure_dir, GRAPH_EGG_TAX_BAR_PNG);
        std::ofstream file_tax_bar(fig_txt_tax_bar, std::ios::out | std::ios::app);
        file_tax_bar << "Taxonomic Scope\tCount" << std::endl;

        stream << "\nTop 10 Taxonomic Scopes Assigned:";
        ct = 1;
        // Sort taxonomy scope
        tax_scope_counter.sort(true);
        for (auto &pair : tax_scope_counter._sorted) {
            if (ct > COUNT_TOP_TAX_SCOPE) break;
            percent = ((fp32)pair.second / tax_scope_counter._ct_total) * 100;
            stream <<
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
        _pGraphingManager->graph(graphingStruct);
    }

    //-------------------------- Gene Ontology --------------------//
    if (ct_total_go_hits > 0) {
        std::string                              fig_txt_bar_go_overall;
        std::string                              fig_png_bar_go_overall;
        std::string                              fig_txt_go_bar;
        std::string                              fig_png_go_bar;

        stream <<
               "\nTotal unique sequences with at least one GO term: " << ct_total_go_hits <<
               "\nTotal unique sequences without GO terms: " << ct_alignments - ct_total_go_hits <<
               "\nTotal GO terms assigned: " << go_combined_map[GO_OVERALL_FLAG]._ct_total;;

        for (uint16 lvl : _go_levels) {
            for (auto &pair : go_combined_map) {
                if (pair.first.empty()) continue;
                // Count maps (biological/molecular/cellular/overall)
                fig_txt_go_bar = PATHS(_figure_dir, pair.first) + std::to_string(lvl)+GRAPH_GO_END_TXT;
                fig_png_go_bar = PATHS(_figure_dir, pair.first) + std::to_string(lvl)+GRAPH_GO_END_PNG;
                std::ofstream file_go_bar(fig_txt_go_bar, std::ios::out | std::ios::app);
                file_go_bar << "Gene Ontology Term\tCount" << std::endl;

                // Sort count maps
                pair.second.sort(true);

                // get total count for each level...change, didn't feel like making another
                uint32 lvl_ct = 0;   // Use for percentages, total terms for each lvl
                ct = 0;              // Use for unique count
                for (auto &pair2 : pair.second._sorted) {
                    if (pair2.first.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                        ct++;
                        lvl_ct += pair2.second;
                    }
                }
                stream << "\nTotal "        << pair.first <<" terms (lvl="          << lvl << "): " << lvl_ct;
                stream << "\nTotal unique " << pair.first <<" terms (lvl="          << lvl << "): " << ct;
                stream << "\nTop 10 "       << pair.first <<" terms assigned (lvl=" << lvl << "): ";

                ct = 1;
                for (auto &pair2 : pair.second._sorted) {
                    if (ct > COUNT_TOP_GO) break;
                    if (pair2.first.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                        percent = ((fp32)pair2.second / lvl_ct) * 100;
                        stream <<
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
                graphingStruct.software_flag = GRAPH_ONTOLOGY_FLAG;
                graphingStruct.graph_type = GRAPH_TOP_BAR_FLAG;
                _pGraphingManager->graph(graphingStruct);
            }
        }
    }

    //--------------------------- KEGG ----------------------------//
    if (ct_total_kegg_hits > 0) {
        stream<<
              "\nTotal unique sequences with at least one pathway (KEGG) assignment: " << ct_total_kegg_hits<<
              "\nTotal unique sequences without pathways (KEGG): " << ct_alignments - ct_total_kegg_hits    <<
              "\nTotal pathways (KEGG) assigned: " << ct_total_kegg_terms;
    }

    out_msg = stream.str();
    _pFileSystem->print_stats(out_msg);
    FS_dprint("Success! EggNOG results parsed");
}

bool ModEggnogDMND::is_executable() {
    std::string test_command;
    TerminalData terminalData;

    test_command = DIAMOND_EXE + " --version";

    terminalData.command = test_command;
    terminalData.print_files = false;

    return TC_execute_cmd(terminalData) == 0;
}

ModEggnogDMND::~ModEggnogDMND() {
    FS_dprint("Killing Object - ModEggnogDMND");

}



std::string ModEggnogDMND::get_output_dmnd_filepath(bool final) {
    std::string filename;

    _blastp ? filename = "blastp" : filename = "blastx";
    filename += "_" + _pUserInput->get_user_transc_basename() + "_eggnog_db";
    if (final) filename += FileSystem::EXT_OUT;
    return PATHS(_mod_out_dir, filename);
}
