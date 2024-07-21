/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2024, Alexander Hart, Dr. Jill Wegrzyn
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
#include "ModEggnog.h"
#include "../ExceptionHandler.h"
#include "../FileSystem.h"
#include "../database/EggnogDatabase.h"

std::vector<ENTAP_HEADERS> ModEggnog::DEFAULT_HEADERS = {
        ENTAP_HEADER_ONT_EGG_SEED_ORTHO,
        ENTAP_HEADER_ONT_EGG_SEED_EVAL,
        ENTAP_HEADER_ONT_EGG_SEED_SCORE,
        ENTAP_HEADER_ONT_EGG_TAX_SCOPE_MAX,
        ENTAP_HEADER_ONT_EGG_MEMBER_OGS,
        ENTAP_HEADER_ONT_EGG_DESC,
        ENTAP_HEADER_ONT_EGG_COG_ABBREVIATION,
        ENTAP_HEADER_ONT_EGG_COG_DESCRIPTION,
        ENTAP_HEADER_ONT_EGG_BIGG,
        ENTAP_HEADER_ONT_EGG_KEGG_KO,
        ENTAP_HEADER_ONT_EGG_KEGG_PATHWAY,
        ENTAP_HEADER_ONT_EGG_KEGG_MODULE,
        ENTAP_HEADER_ONT_EGG_KEGG_REACTION,
        ENTAP_HEADER_ONT_EGG_KEGG_RCLASS,
        ENTAP_HEADER_ONT_EGG_BRITE,
        ENTAP_HEADER_ONT_EGG_GO_BIO,
        ENTAP_HEADER_ONT_EGG_GO_CELL,
        ENTAP_HEADER_ONT_EGG_GO_MOLE,
        ENTAP_HEADER_ONT_EGG_PROTEIN
};


ModEggnog::ModEggnog(std::string &ont_out, std::string &in_hits, EntapDataPtrs &entap_data)
        : AbstractOntology(in_hits, ont_out, entap_data, "EggNOG", DEFAULT_HEADERS) {

    mExePath = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_MAPPER_EXE);
    mEggnogMapDMNDPath = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_MAPPER_DMND_DB);
    mEggnogMapDataDir = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_MAPPER_DATA_DIR);

    mEggnogMapAnnotationsOutputPath = PATHS(mModOutDir, get_output_tag()+EGG_OUTPUT_ANNOT_APPEND);
    mEggnogMapHitsOutputPath = PATHS(mModOutDir, get_output_tag() + EGG_OUTPUT_HITS_APPEND);
    mEggnogMapSeedOrthoOutputPath = PATHS(mModOutDir, get_output_tag() + EGG_OUTPUT_SEED_ORTHO_APPEND);
    mEggnogMapperState = EGGNOG_MAPPER_NOT_STARTED;
    mSoftwareFlag = ONT_EGGNOG_MAPPER;
}

bool ModEggnog::is_executable(std::string &exe) {
    TerminalData terminalData;

    terminalData.command = exe + " --version";
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;
    return TC_execute_cmd(terminalData) == 0;
}


EntapModule::ModVerifyData ModEggnog::verify_files() {
    ModVerifyData modVerifyData;
    modVerifyData.files_exist = false;
    uint16 file_status;

    FS_dprint("Overwrite was unselected, verifying output files...");


    // Eggnog-mapper generates several files, we are mainly concerned with the annotations and seed orthologs
    //  The 'TAG' can be set when running EggNOG mapper to have a specific output tag in the filename
    //  1. TAG.emapper.annotations
    //  2. TAG.emapper.hits
    //  3. TAG.emapper.seed_orthologs

    file_status = mpFileSystem->get_file_status(mEggnogMapAnnotationsOutputPath);
    if (file_status == 0) {
        FS_dprint("Valid annotations file found at: " + mEggnogMapAnnotationsOutputPath);
        mEggnogMapperState = EGGNOG_MAPPER_ANNOTATIONS_COMPLETE;
        modVerifyData.files_exist = true;
    } else {
        // Annotations file not found, look for seed orthologs
        file_status = mpFileSystem->get_file_status(mEggnogMapSeedOrthoOutputPath);
        if (file_status == 0) {
            FS_dprint("Valid seed ortholog file found at: " + mEggnogMapSeedOrthoOutputPath);
            mpFileSystem->delete_file(mEggnogMapAnnotationsOutputPath);
            mEggnogMapperState = EGGNOG_MAPPER_DIAMOND_COMPLETE;
        } else {
            mpFileSystem->delete_file(mEggnogMapAnnotationsOutputPath);
            mpFileSystem->delete_file(mEggnogMapSeedOrthoOutputPath);
            FS_dprint("No Eggnog-mapper output files found, running...");
        }
    }
    return modVerifyData;
}


/**
 * ======================================================================
 * Function void ModEggnog::execute()
 *
 * Description          - Main execution routine
 *                      - Generates command for pstreams and runs eggnog
 *                        against DIAMOND hits and DIAMOND no-hits
 *
 * Notes                - None
 *
 *
 * @return              - None
 * ======================================================================
 */
void ModEggnog::execute() {
    std::string                        std_out;
    std::string                        cmd;
    command_map_t                      tc_command_map;
    TerminalData                       terminalData;

    // Ensure we do not have an annotations file, should never occur (checked in other function)
    if (mpFileSystem->file_exists(mEggnogMapAnnotationsOutputPath)) return;
    TC_print(TC_PRINT_COUT, "Running EggNOG-mapper analysis...");

    tc_command_map = {
            {"-m", "diamond"},
            {"-o", get_output_tag()},
            {"--sensmode", "more-sensitive"},
            {"-i", mInputTranscriptome},
            {"--dmnd_db", mEggnogMapDMNDPath},
            {"--data_dir", mEggnogMapDataDir},
            {"--output_dir", mModOutDir},
            {"--cpu", std::to_string(mpUserInput->get_supported_threads())},
            {"--no_file_comments",""}
    };

    // If we are in runN/blastx mode, run eggNOG as such
    if (!mBlastp) {
        tc_command_map.emplace("--itype", "CDS");
    }

    // Does user want to specify --dbmem flag for EggNOG
    if (mpUserInput->has_input(INPUT_FLAG_EGG_MAPPER_DBMEM))
    {
        tc_command_map.emplace("--dbmem", "");
    }

    switch (mEggnogMapperState) {
        case EGGNOG_MAPPER_NOT_STARTED:
            FS_dprint("Running EggNOG-mapper from beginning...");
            break;

        case EGGNOG_MAPPER_DIAMOND_COMPLETE:
            FS_dprint("Running EggNOG-mapper resuming after DIAMOND completion...");
            tc_command_map.emplace("--resume","");
            // Other commands can remain the same when we are resuming
            break;

        case EGGNOG_MAPPER_ANNOTATIONS_COMPLETE:
        default:
            FS_dprint("ERROR unhandled state in ModEggnog execute");
            return; // WARNING!!!!! returning from invalid state
            break;
    }

    cmd = TC_generate_command(tc_command_map, mExePath);
    std_out = PATHS(mModOutDir, get_output_tag()) + "_" + FileSystem::EXT_STD;
    terminalData.command = cmd;
    terminalData.suppress_std_err = false;
    terminalData.base_std_path = std_out;
    terminalData.print_files = true;

    // Execute Eggnog and WAIT for completion
    if (TC_execute_cmd(terminalData) != TC_EXIT_SUCCESS) {
        // Clear all output files for now, it may be possible for EggNOG to complete the hits stage and fail
        //  annotations stage so we can do an improvement here. Right now all or none though. Errcode from eggnog looks like always 1
        mpFileSystem->delete_file(mEggnogMapAnnotationsOutputPath);
        mpFileSystem->delete_file(mEggnogMapSeedOrthoOutputPath);
        mpFileSystem->delete_file(mEggnogMapHitsOutputPath);

        throw ExceptionHandler("Error in running EggNOG-mapper\nEggNOG Error:\n" + terminalData.err_stream,
                               ERR_ENTAP_RUN_EGGNOG);
    }

    TC_print(TC_PRINT_COUT, "Success!");
    FS_dprint("Success! EggNOG execution complete");
}


/**
 * ======================================================================
 * Function void ModEggnog::parse()
 *
 * Description          - Analyzes eggnog output, parsing and graphing results
 *
 * Notes                - None
 *
 *
 * @return              - None
 * ======================================================================
 */
void ModEggnog::parse() {

    FS_dprint("Beginning to parse eggnog mapper results...");
    TC_print(TC_PRINT_COUT, "Parsing EggNOG-mapper Analysis...");

    std::string                              out_msg;
    uint32                                   count_total_go_hits=0;
    uint32                                   count_TOTAL_hits=0;         // All ortho matches
    uint32                                   count_total_kegg_hits=0;
    uint32                                   count_no_hits=0;            // Unannotated OGs
    uint32                                   ct = 0;
    fp32                                     percent;
    Compair<std::string>                     tax_scope_counter;
    std::vector<ENTAP_HEADERS>               output_headers;
    std::unordered_map<std::string,Compair<GoEntry>>  go_combined_map;     // Just for convenience
    QuerySequence                           *querySequence;
    QuerySequence::EggnogResults            EggnogResults;
    EggnogDatabase                          *pEggnogDatabase;
    std::stringstream stats_stream;

    FS_dprint("Eggnog file located at " + mEggnogMapAnnotationsOutputPath + " being filtered");
    if (!mpFileSystem->file_exists(mEggnogMapAnnotationsOutputPath)) {
        throw ExceptionHandler("EggNOG file not found at: " + mEggnogMapAnnotationsOutputPath, ERR_ENTAP_PARSE_EGGNOG);
    }

    // setup headers for printing
    output_headers = DEFAULT_HEADERS;
    output_headers.insert(output_headers.begin(), ENTAP_HEADER_QUERY);

    // Get EggNOG database, only used to map COG categories currently
    pEggnogDatabase = new EggnogDatabase(mpFileSystem, mpEntapDatabase, mpQueryData);

    // Begin to read through TSV file, these are all the headers in a  default eggnog-mapper run
    std::string qseqid, seed_ortho, seed_score, eggnog_ogs, max_annot_tax_level, cog_category, description,
            preferred_name, gene_ontology_terms, ec_value, kegg_ko, kegg_pathway, kegg_mod, kegg_reaction, kegg_rclass,
            brite, kegg_tc, cazy, bigg_reaction, pfams;
    fp64 seed_e;
    io::CSVReader<EGGNOG_COL_NUM, io::trim_chars<' ','-'>, io::no_quote_escape<'\t'>> in(mEggnogMapAnnotationsOutputPath);
    in.next_line(); // Skip header line
    while (in.read_row(qseqid, seed_ortho, seed_e, seed_score, eggnog_ogs, max_annot_tax_level, cog_category, description,
                       preferred_name, gene_ontology_terms, ec_value, kegg_ko, kegg_pathway, kegg_mod, kegg_reaction, kegg_rclass,
                       brite, kegg_tc, cazy, bigg_reaction, pfams)) {

        // Check if the query matches one of our original transcriptome sequences
        querySequence = mpQueryData->get_sequence(qseqid);
        if (querySequence == nullptr) {
            delete pEggnogDatabase;
            throw ExceptionHandler("Unable to find sequence " + qseqid + " in input transcriptome",
                                   ERR_ENTAP_PARSE_EGGNOG);
        }

        count_TOTAL_hits++;

        // Populate data from EggNOG-mapper run
        EggnogResults = {};
        EggnogResults.seed_ortholog = seed_ortho;
        EggnogResults.seed_evalue = float_to_sci(seed_e, 2);
        EggnogResults.seed_eval_raw = seed_e;
        EggnogResults.seed_score = seed_score;
        EggnogResults.member_ogs = eggnog_ogs;
        EggnogResults.tax_scope_lvl_max = max_annot_tax_level;
        EggnogResults.cog_category = cog_category;
        // Ensure COG Category abbreviation is not empty or NULL
        if ((!cog_category.empty()) && (cog_category != EGGNOG_NULL_CHARACTER)) {
            // Format from EggNOG is "C" or "CA", can be multiple
            std::string cog_descriptions;
            for (char abbrev : cog_category) {
                cog_descriptions += pEggnogDatabase->get_cog_category_description(abbrev) + ";";
            }
            cog_descriptions.pop_back();
            EggnogResults.cog_category_description = cog_descriptions;
        }
        EggnogResults.description = description;
        EggnogResults.pname = preferred_name;
        EggnogResults.parsed_go = mpEntapDatabase->format_go_delim(gene_ontology_terms,',');
        EggnogResults.ec_value = ec_value;
        EggnogResults.kegg_ko = kegg_ko;
        EggnogResults.kegg_pathway = kegg_pathway;
        EggnogResults.kegg_module = kegg_mod;
        EggnogResults.kegg_reaction = kegg_reaction;
        EggnogResults.kegg_rclass = kegg_rclass;
        EggnogResults.brite = brite;
        EggnogResults.kegg_tc = kegg_tc;
        EggnogResults.cazy = cazy;
        EggnogResults.bigg = bigg_reaction;
        EggnogResults.protein_domains = pfams;

        querySequence->add_alignment(GENE_ONTOLOGY, mSoftwareFlag, EggnogResults, mEggnogMapDMNDPath);

        // Count GO results
        if (!EggnogResults.parsed_go.empty()) {
            count_total_go_hits++;
            for (auto goEntry : EggnogResults.parsed_go) {
                go_combined_map[goEntry.category].add_value(goEntry);
                go_combined_map[GO_OVERALL_FLAG].add_value(goEntry);
            }
        }

        // Count KEGG stats
        if (querySequence->QUERY_FLAG_GET(QuerySequence::QUERY_FAMILY_ONE_KEGG)) {
            count_total_kegg_hits++;
        }

        // Compile Taxonomic Orthogroup stats
        if (!EggnogResults.tax_scope_lvl_max.empty()) {
            // Count number of unique taxonomic groups
            tax_scope_counter.add_value(EggnogResults.tax_scope_lvl_max);
        }

    } // END WHILE file reading
    delete pEggnogDatabase;

    FS_dprint("Success! Printing output files");
    // Initialize and print output files, inefficient redo
    uint32 output_annotated_query_flags=0;
    output_annotated_query_flags |= QuerySequence::QUERY_EGGNOG_HIT;
    output_annotated_query_flags |= QuerySequence::QUERY_FRAME_KEPT;
    output_annotated_query_flags |= QuerySequence::QUERY_EXPRESSION_KEPT;
    uint32 output_unannotated_query_flags = 0;
    output_unannotated_query_flags |= QuerySequence::QUERY_FRAME_KEPT;
    output_unannotated_query_flags |= QuerySequence::QUERY_EXPRESSION_KEPT;
    if (mpQueryData->is_protein_data()) {
        std::string out_no_hits_faa = PATHS(mProcDir, EGG_MAPPER_PREFIX + EGG_OUT_UNANNOTATED + FileSystem::EXT_FAA);
        mpQueryData->print_transcriptome(output_unannotated_query_flags, out_no_hits_faa, QueryData::SEQUENCE_AMINO_ACID);
        std::string out_hits_faa = PATHS(mProcDir, EGG_MAPPER_PREFIX + EGG_OUT_ANNOTATED + FileSystem::EXT_FAA);
        mpQueryData->print_transcriptome(output_annotated_query_flags, out_hits_faa, QueryData::SEQUENCE_AMINO_ACID);
    }

    if (mpQueryData->is_nucleotide_data()) {
        std::string out_no_hits_fnn = PATHS(mProcDir, EGG_MAPPER_PREFIX + EGG_OUT_UNANNOTATED + FileSystem::EXT_FNN);
        std::string out_hits_fnn = PATHS(mProcDir, EGG_MAPPER_PREFIX + EGG_OUT_ANNOTATED + FileSystem::EXT_FNN);
        mpQueryData->print_transcriptome(output_unannotated_query_flags, out_no_hits_fnn, QueryData::SEQUENCE_NUCLEOTIDE);
        mpQueryData->print_transcriptome(output_annotated_query_flags, out_hits_fnn, QueryData::SEQUENCE_NUCLEOTIDE);
    }

    count_no_hits = mpQueryData->getMTotalSequences() - count_TOTAL_hits;

    FS_dprint("Success! Computing overall statistics...");

    // Begin stats stream
    mpFileSystem->format_stat_stream(stats_stream, "Gene Family - Gene Ontology and Pathway - EggNOG");
    // Make sure we have any alignmnts against EggNOG database
    if (count_TOTAL_hits > 0) {
        stats_stream <<
               "Statistics for overall Eggnog results: "               <<
               "\nTotal unique sequences with family assignment: "     << count_TOTAL_hits <<
               "\nTotal unique sequences without family assignment: "  << count_no_hits;

        //--------------------- Top Ten Taxonomic Scopes --------------//
        if (!tax_scope_counter.empty()) {
            stats_stream << "\nTop " << std::to_string(COUNT_TOP_TAX_SCOPE) << " Taxonomic Scopes Assigned:";
            ct = 1;
            // Sort taxonomy scope
            tax_scope_counter.sort(true);
            for (auto &pair : tax_scope_counter._sorted) {
                if (ct > COUNT_TOP_TAX_SCOPE) break;
                percent = ((fp32)pair.second / (fp32)mpQueryData->getMTotalKeptSequences()) * 100;
                stats_stream <<
                       "\n\t" << ct << ")" << pair.first << ": " << pair.second <<
                       "(" << percent << "% of total retained sequences)";
                ct++;
            }
        }
        //-------------------------------------------------------------//

        //-------------------------- Gene Ontology --------------------//
        if (count_total_go_hits > 0) {
            std::string                              fig_txt_bar_go_overall;
            std::string                              fig_png_bar_go_overall;
            std::string                              fig_txt_go_bar;
            std::string                              fig_png_go_bar;

            stats_stream <<
                   "\nTotal unique sequences with at least one GO term: " << count_total_go_hits <<
                   "\nTotal unique sequences without GO terms: " << count_TOTAL_hits - count_total_go_hits <<
                   "\nTotal GO terms assigned: " << go_combined_map[GO_OVERALL_FLAG]._ct_total;

            for (auto &pair : go_combined_map) {
                if (pair.first.empty() || pair.second.empty()) continue;
                // Sort count maps
                pair.second.sort(true);

                // get total count for each category
                uint32 lvl_ct = 0;   // Use for percentages, total terms for each lvl
                ct = 0;              // Use for unique count
                // pair2 = unique go entry : number of that GoEntry
                for (auto &pair2 : pair.second._sorted) {
                    ct++;
                    lvl_ct += pair2.second;
                }
                stats_stream << "\nTotal "        << pair.first <<" terms: " << lvl_ct;
                stats_stream << "\nTotal unique " << pair.first <<" terms: " << ct;
                stats_stream << "\nTop " << COUNT_TOP_GO << " " << pair.first <<" terms assigned: ";

                // Get the TOP x go terms and print out based on occurance
                ct = 1;
                for (auto &pair2 : pair.second._sorted) {
                    if (ct > COUNT_TOP_GO) break;
                    percent = ((fp32)pair2.second / lvl_ct) * 100;
                    stats_stream <<
                           "\n\t" << ct << ")" << pair2.first.go_id << ": " << pair2.second <<
                           "(" << percent << "% of total GO terms)";
                    ct++;
                }
            }
        }
        //-------------------------------------------------------------//

        //--------------------------- KEGG ----------------------------//
        if (count_total_kegg_hits > 0) {
            stats_stream <<
                  "\nTotal unique sequences with at least one KEGG assignment: " << count_total_kegg_hits<<
                  "\nTotal unique sequences without KEGG assignment: " << count_TOTAL_hits - count_total_kegg_hits;
        }
        //-------------------------------------------------------------//

    } else {
        FS_dprint("WARNING: NO alignments against EggNOG!");
        stats_stream << "Warning: No alignments against EggNOG database" << std::endl;
    }

    out_msg = stats_stream.str();
    mpFileSystem->print_stats(out_msg);
    FS_dprint("Success! EggNOG results parsed");
}

ModEggnog::~ModEggnog() {
    FS_dprint("Killing object - ModEggnog");
}

std::string ModEggnog::get_output_tag() {
    std::string ret_tag;
    mBlastp ? ret_tag = "blastp_" : ret_tag = "blastx_";
    ret_tag += mpUserInput->get_user_transc_basename();
    return ret_tag;
}

bool ModEggnog::set_version() {
    return false;
}
