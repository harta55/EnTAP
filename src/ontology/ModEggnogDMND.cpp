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

#include <csv.h>
#include "ModEggnogDMND.h"
#include "../database/EggnogDatabase.h"
#include "../TerminalCommands.h"
#include "../QueryAlignment.h"
#include "../QueryData.h"
#include "../GraphingManager.h"
#include "../ExceptionHandler.h"


std::vector<ENTAP_HEADERS> ModEggnogDMND::DEFAULT_HEADERS = {
    ENTAP_HEADER_ONT_EGG_SEED_ORTHO,
    ENTAP_HEADER_ONT_EGG_SEED_EVAL,
    ENTAP_HEADER_ONT_EGG_SEED_SCORE,
    ENTAP_HEADER_ONT_EGG_PRED_GENE,
    ENTAP_HEADER_ONT_EGG_TAX_SCOPE_READABLE,
    ENTAP_HEADER_ONT_EGG_TAX_SCOPE_MAX,
    ENTAP_HEADER_ONT_EGG_MEMBER_OGS,
    ENTAP_HEADER_ONT_EGG_DESC,
    ENTAP_HEADER_ONT_EGG_KEGG,
    ENTAP_HEADER_ONT_EGG_BIGG,
    ENTAP_HEADER_ONT_EGG_GO_BIO,
    ENTAP_HEADER_ONT_EGG_GO_CELL,
    ENTAP_HEADER_ONT_EGG_GO_MOLE,
    ENTAP_HEADER_ONT_EGG_PROTEIN
};

ModEggnogDMND::ModEggnogDMND(std::string &ont_out, std::string &in_hits,
                             EntapDataPtrs &entap_data)
        : AbstractOntology(in_hits, ont_out, entap_data, "EggNOG_DMND", DEFAULT_HEADERS) {
    FS_dprint("Spawn Object - ModEggnogDMND");

    mSoftwareFlag = ONT_EGGNOG_DMND;
    mEggnogDbDiamond = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_DMND_DB);
    mEggnogDbSQL     = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_SQL_DB);
    mExePath         = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_DIAMOND_EXE);
}

EntapModule::ModVerifyData ModEggnogDMND::verify_files() {
    ModVerifyData modVerifyData;
    modVerifyData.files_exist = false;
    uint16 file_status = 0;

    FS_dprint("Overwrite was unselected, verifying output files...");
    mOutHIts = get_output_dmnd_filepath(true);
    file_status = mpFileSystem->get_file_status(mOutHIts);

    if (file_status != 0) {
        FS_dprint(mpFileSystem->print_file_status(file_status,mOutHIts));
        FS_dprint("Errors in opening file, continuing with execution...");
        modVerifyData.files_exist = false;
    } else {
        modVerifyData.files_exist = true;
    }
    modVerifyData.output_paths = vect_str_t{mOutHIts};
    return modVerifyData;
}

void ModEggnogDMND::execute() {
    std::string                        std_out;
    std::string                        cmd;
    std::string                        blast;
    TerminalData                       terminalData;

    FS_dprint("Running EggNOG against Diamond database...");

    // Ensure both input path and EggNOG DMND database exist before continuing
    if (!mpFileSystem->file_exists(mEggnogDbDiamond)) {
        throw ExceptionHandler("EggNOG DIAMOND database not found at: " + mEggnogDbDiamond,
                               ERR_ENTAP_EGGNOG_FILES);
    }
    if (!mpFileSystem->file_exists(mInputTranscriptome)) {
        throw ExceptionHandler("Input transcriptome not found at: " + mInputTranscriptome, ERR_ENTAP_EGGNOG_FILES);
    }

    // Generate paths for DIAMOND run (out_hits set previously)
    std_out = get_output_dmnd_filepath(false) + "_" + FileSystem::EXT_STD;
    mBlastp ? blast = "blastp" : blast = "blastx";

    //Run DIAMOND
    cmd =
            mExePath + " " +
            blast +
            " -d " + mEggnogDbDiamond +
            " --top 1"             +
            " --more-sensitive"    +
            " -q "                 + mInputTranscriptome   +
            " -o "                 + mOutHIts  +
            " -p "                 + std::to_string(mThreads) +
            " -f " + "6 qseqid sseqid pident length mismatch gapopen "
                    "qstart qend sstart send evalue bitscore qcovhsp stitle";

    terminalData.command        = cmd;
    terminalData.print_files    = true;
    terminalData.suppress_std_err = false;
    terminalData.base_std_path  = std_out;

    if (TC_execute_cmd(terminalData) != 0) {
        // Error in run
        mpFileSystem->delete_file(mOutHIts);
        FS_dprint("DIAMOND STD OUT:\n" + terminalData.out_stream);
        throw ExceptionHandler("Error in running DIAMOND against EggNOG database at: " +
                               mEggnogDbDiamond + "\nDIAMOND Error:\n" + terminalData.err_stream, ERR_ENTAP_RUN_EGGNOG_DMND);
    }
}

void ModEggnogDMND::parse() {
    uint16         file_status=0;
    uint64         sequence_ct=0;   // dprintf sequence count
    std::stringstream stats_stream;

    FS_dprint("Parsing EggNOG DMND file located at: " + mOutHIts);

    // Ensure file is valid
    file_status = mpFileSystem->get_file_status(mOutHIts);
    if (file_status != 0) {
        throw ExceptionHandler(mpFileSystem->print_file_status(file_status,mOutHIts),
                               ERR_ENTAP_PARSE_EGGNOG_DMND);
    }

    // File valid, continue
    FS_dprint("Beginning to parse EggNOG results...");
    mpFileSystem->format_stat_stream(stats_stream, "Gene Family - Gene Ontology and Pathway - EggNOG");

#ifdef USE_FAST_CSV
    // ------------------ Read from DIAMOND output ---------------------- //
    std::string qseqid;
    std::string sseqid, stitle, database_name,pident, bitscore,
            length, mismatch, gapopen, qstart, qend, sstart, send, coverage;
    fp64 evalue;
    QuerySequence::EggnogResults eggnogResults;
    QuerySequence *querySequence;
    // ----------------------------------------------------------------- //
    // Begin using CSVReader lib to parse data
    try {
        io::CSVReader<DMND_COL_NUMBER, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(mOutHIts);
        while (in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                           qstart, qend, sstart, send, evalue, bitscore, coverage,stitle)) {
            // Currently throwing away most DIAMOND results

            // Print progress to debug
            if (++sequence_ct % STATUS_UPDATE_HITS == 0) {
                FS_dprint("Alignments parsed: " + std::to_string(sequence_ct));
            }

            // Ensure we recognize the query sequence before continuing
            querySequence = mpQueryData->get_sequence(qseqid);
            if (querySequence == nullptr) {
                throw ExceptionHandler("Unable to find sequence " + qseqid + " in input transcriptome",
                                       ERR_ENTAP_PARSE_EGGNOG_DMND);
            }

            // Populate seed data from diamond
            eggnogResults = {};
            eggnogResults.seed_eval_raw = evalue;
            eggnogResults.seed_evalue = float_to_sci(evalue, 2);
            eggnogResults.seed_score  = bitscore;
            eggnogResults.seed_coverage = coverage;
            eggnogResults.seed_ortholog = sseqid;

            // WARNING!!! SQL lookups are done in "calculate_stats" below to save execution time
            //      (only best hits are looked up) headers are populated then!
            querySequence->add_alignment(GENE_ONTOLOGY, mSoftwareFlag, eggnogResults, mEggnogDbDiamond);

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
    FS_dprint("Success! Calculating statistics and accessing EggNOG database...");

    QuerySequence::EggnogResults                          *eggnog_results;
    EggnogDmndAlignment                                  *best_hit;
    Compair<std::string>                                  tax_scope_counter;
    std::unordered_map<std::string,Compair<GoEntry>>  go_combined_map;
    EggnogDatabase                                       *eggnogDatabase;
    std::vector<ENTAP_HEADERS>                            output_headers;
    GraphingManager::GraphingData                         graphing_data_temp;

    uint64         ct_alignments=0;
    uint64         ct_no_alignment=0;
    uint64         ct_total_go_hits=0;      // Sequences that had at least one go
    uint64         ct_total_kegg_hits=0;    // Sequences that had at least one kegg
    uint64         ct_total_kegg_terms=0;
    uint32         ct = 0;
    fp32           percent;

    std::string    out_msg;

    // setup headers for printing
    output_headers = DEFAULT_HEADERS;
    output_headers.insert(output_headers.begin(), ENTAP_HEADER_QUERY);

    // Generate EggNOG database
    eggnogDatabase = new EggnogDatabase(mpFileSystem, mpEntapDatabase, mpQueryData);
    if (eggnogDatabase->open_sql(mEggnogDbSQL) != EggnogDatabase::ERR_EGG_OK) {
        delete eggnogDatabase;
        throw ExceptionHandler("Unable to open EggNOG SQL Database", ERR_ENTAP_PARSE_EGGNOG_DMND);
    }

    // Output files
    std::string out_no_hits_base = PATHS(mProcDir, FILENAME_OUT_UNANNOTATED);
    std::string out_hits_base    = PATHS(mProcDir, FILENAME_OUT_ANNOTATED);

    mpQueryData->start_alignment_files(out_no_hits_base, output_headers, mGoLevels, mAlignmentFileTypes);
    mpQueryData->start_alignment_files(out_hits_base, output_headers, mGoLevels, mAlignmentFileTypes);

    // Parse through all query sequences
    for (auto &pair : *mpQueryData->get_sequences_ptr()) {
        // Check if each sequence is an eggnog alignment
        if (pair.second->hit_database(GENE_ONTOLOGY, mSoftwareFlag, mEggnogDbDiamond)) {
            // Yes, hit EggNOG database
            ct_alignments++;

            best_hit = pair.second->get_best_hit_alignment<EggnogDmndAlignment>
                    (GENE_ONTOLOGY, mSoftwareFlag, mEggnogDbDiamond);

            if (best_hit == nullptr) {
                FS_dprint("ERROR nullptr pulled ModEggnogDMND: " + pair.first);
                continue;
            }
            eggnog_results = best_hit->get_results();
            eggnogDatabase->get_eggnog_entry(eggnog_results);
            best_hit->refresh_headers();

            mpQueryData->add_alignment_data(out_hits_base, pair.second, nullptr);

            //  Analyze Gene Ontology Stats
            if (!eggnog_results->parsed_go.empty()) {
                ct_total_go_hits++;
                for (auto &go_entry: eggnog_results->parsed_go) {
                    // pair - first: GO category, second; vector of terms
                    go_combined_map[go_entry.category].add_value(go_entry);
                    go_combined_map[GO_OVERALL_FLAG].add_value(go_entry);
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
            mpQueryData->add_alignment_data(out_no_hits_base, pair.second, nullptr);
        }
    } // END FOR LOOP

    // Close files
    mpQueryData->end_alignment_files(out_hits_base);
    mpQueryData->end_alignment_files(out_no_hits_base);
    delete eggnogDatabase;

    FS_dprint("EggNOG database closed, printing stats...");

    // Begin to print stats / files
    stream <<
       "Statistics for overall Eggnog results: "               <<
       "\nTotal unique sequences with family assignment: "     << ct_alignments <<
       "\nTotal unique sequences without family assignment: "  << ct_no_alignment;

    // Make sure we have hits before doing anything

    //--------------------- Top Ten Taxonomic Scopes --------------//
    if (!tax_scope_counter.empty()) {
        // Setup graphing files
        graphing_data_temp = GraphingManager::GraphingData();
        graphing_data_temp.x_axis_label = "Taxonomic Scope";
        graphing_data_temp.y_axis_label = "Count";
        graphing_data_temp.text_file_path = PATHS(mFigureDir, GRAPH_EGG_TAX_BAR_TXT);
        graphing_data_temp.fig_out_path   = PATHS(mFigureDir, GRAPH_EGG_TAX_BAR_PNG);
        graphing_data_temp.graph_title    = GRAPH_EGG_TAX_BAR_TITLE;
        graphing_data_temp.graph_type     = GraphingManager::ENT_GRAPH_BAR_HORIZONTAL;

        mpGraphingManager->initialize_graph_data(graphing_data_temp);

        stream << "\nTop " << std::to_string(COUNT_TOP_TAX_SCOPE) << " Taxonomic Scopes Assigned:";
        ct = 1;
        // Sort taxonomy scope
        tax_scope_counter.sort(true);
        for (auto &pair : tax_scope_counter._sorted) {
            if (ct > COUNT_TOP_TAX_SCOPE) break;
            percent = ((fp32)pair.second / tax_scope_counter._ct_total) * 100;
            stream <<
               "\n\t" << ct << ")" << pair.first << ": " << pair.second <<
               "(" << percent << "%)";
            mpGraphingManager->add_datapoint(graphing_data_temp.text_file_path, {pair.first, std::to_string(pair.second)});
            ct++;
        }
        mpGraphingManager->graph_data(graphing_data_temp.text_file_path);
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

        for (uint16 lvl : mGoLevels) {
            for (auto &pair : go_combined_map) {
                if (pair.first.empty() || pair.second.empty()) continue;
                // Count maps (biological/molecular/cellular/overall)
                graphing_data_temp = GraphingManager::GraphingData();
                graphing_data_temp.x_axis_label = "Gene Ontology Term";
                graphing_data_temp.y_axis_label = "Count";
                graphing_data_temp.text_file_path = PATHS(mFigureDir, pair.first) + std::to_string(lvl)+GRAPH_GO_END_TXT;
                graphing_data_temp.fig_out_path   = PATHS(mFigureDir, pair.first) + std::to_string(lvl)+GRAPH_GO_END_PNG;

                if (pair.first == GO_BIOLOGICAL_FLAG) graphing_data_temp.graph_title = GRAPH_GO_BAR_BIO_TITLE + "_Level:_"+std::to_string(lvl);
                if (pair.first == GO_CELLULAR_FLAG) graphing_data_temp.graph_title = GRAPH_GO_BAR_CELL_TITLE+ "_Level:_"+std::to_string(lvl);
                if (pair.first == GO_MOLECULAR_FLAG) graphing_data_temp.graph_title = GRAPH_GO_BAR_MOLE_TITLE+ "_Level:_"+std::to_string(lvl);
                if (pair.first == GO_OVERALL_FLAG) graphing_data_temp.graph_title = GRAPH_GO_BAR_ALL_TITLE+ "_Level:_"+std::to_string(lvl);
                graphing_data_temp.graph_type = GraphingManager::ENT_GRAPH_BAR_HORIZONTAL;

                mpGraphingManager->initialize_graph_data(graphing_data_temp);

                // Sort count maps
                pair.second.sort(true);

                // get total count for each level
                uint32 lvl_ct = 0;   // Use for percentages, total terms for each lvl
                ct = 0;              // Use for unique count
                // pair2 = unique go entry : number of that GoEntry
                for (auto &pair2 : pair.second._sorted) {
                    if (pair2.first.level_int >= lvl || lvl == 0) {
                        ct++;
                        lvl_ct += pair2.second;
                    }
                }
                stream << "\nTotal "        << pair.first <<" terms (lvl>="          << lvl << "): " << lvl_ct;
                stream << "\nTotal unique " << pair.first <<" terms (lvl>="          << lvl << "): " << ct;
                stream << "\nTop " << COUNT_TOP_GO << " " << pair.first <<" terms assigned (lvl>=" << lvl << "): ";

                // Get the TOP x go terms and print out based on occurance
                ct = 1;
                for (auto &pair2 : pair.second._sorted) {
                    if (ct > COUNT_TOP_GO) break;
                    if (pair2.first.level_int >= lvl || lvl == 0) {
                        percent = ((fp32)pair2.second / lvl_ct) * 100;
                        stream <<
                               "\n\t" << ct << ")" << pair2.first.go_id << ": " << pair2.second <<
                               "(" << percent << "%)";
                        mpGraphingManager->add_datapoint(graphing_data_temp.text_file_path, {pair2.first.go_id, std::to_string(pair2.second)});
                        ct++;
                    }
                }
                mpGraphingManager->graph_data(graphing_data_temp.text_file_path);
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
    mpFileSystem->print_stats(out_msg);
    FS_dprint("Success! EggNOG results parsed");
}

bool ModEggnogDMND::is_executable(std::string &exe) {
    std::string test_command;
    TerminalData terminalData;

    test_command = exe + " --version";

    terminalData.command = test_command;
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;

    return TC_execute_cmd(terminalData) == 0;
}

ModEggnogDMND::~ModEggnogDMND() {
    FS_dprint("Killing Object - ModEggnogDMND");
}

std::string ModEggnogDMND::get_output_dmnd_filepath(bool final) {
    std::string filename;

    mBlastp ? filename = "blastp" : filename = "blastx";
    filename += "_" + mpUserInput->get_user_transc_basename() + "_eggnog_proteins";
    if (final) filename += FileSystem::EXT_OUT;
    return PATHS(mModOutDir, filename);
}

bool ModEggnogDMND::set_version() {
    return false;

}
