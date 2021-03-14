#ifdef EGGNOG_MAPPER

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


#include <boost/archive/binary_iarchive.hpp>
#include <csv.h>
#include "ModEggnog.h"
#include "../ExceptionHandler.h"
#include "../FileSystem.h"

std::string ModEggnog::EGG_EMAPPER_EXE = "";

/**
 * ======================================================================
 * Function std::pair<bool, std::string> ModEggnog::verify_files()
 *
 * Description          - Checks whether eggnog has already been ran
 *                      - Will not be called if 'overwrite' is selected by user
 *
 * Notes                - None
 *
 *
 * @return              - Pair (first: files found, second: unused)
 * ======================================================================
 */
std::pair<bool, std::string> ModEggnog::verify_files() {
    std::string                        annotation_base_flag;

    annotation_base_flag = PATHS(_egg_out_dir, EGG_ANNOT_RESULTS);
    mOutHIts            = annotation_base_flag  + EGG_ANNOT_APPEND;

    FS_dprint("Overwrite was unselected, verifying output files...");
    if (pFileSystem->file_exists(mOutHIts)) {
        FS_dprint("File located at: " + mOutHIts + " found, skipping EggNOG execution");
        return std::make_pair(true, "");
    } else {
        FS_dprint("File located at: " + mOutHIts + " NOT found, continuing execution");
        return std::make_pair(false, "");
    }
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

    FS_dprint("Beginning to parse eggnog results...");

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
    std::string                              fig_txt_go_bar;
    std::string                              fig_png_go_bar;
    uint32                                   count_total_go_hits=0;
    uint32                                   count_no_go=0;
    uint32                                   count_no_kegg=0;
    uint32                                   count_TOTAL_hits=0;         // All ortho matches
    uint32                                   count_total_kegg_terms=0;
    uint32                                   count_total_kegg_hits=0;
    uint32                                   count_no_hits=0;            // Unannotated OGs
    uint32                                   ct = 0;
    fp32                                     percent;
    SQLDatabaseHelper                        EGGNOG_DATABASE;
    Compair<std::string>                     tax_scope_counter;
    std::unordered_map<std::string,Compair<std::string>>  go_combined_map;     // Just for convenience
    GraphingData                             graphingStruct;

    QuerySequence::EggnogResults            EggnogResults;

    ss<<std::fixed<<std::setprecision(2);

    if (!EGGNOG_DATABASE.open(mEggnogDbPath))
        throw ExceptionHandler("Unable to open EggNOG database",ERR_ENTAP_PARSE_EGGNOG);

    path = mOutHIts;

    FS_dprint("Eggnog file located at " + path + " being filtered");
    if (!pFileSystem->file_exists(path)) {
        FS_dprint("File not found, at: " + path);
        throw ExceptionHandler("EggNOG file not found at: " + path, ERR_ENTAP_PARSE_EGGNOG);
    }

    // Format EggNOG output - remove comments (the EggNOG mapper command doesn't work correctly), make sure all columns are equal
    path = eggnog_format(path);

    // Begin to read through TSV file
    std::string qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
            best_og, cog_cat, eggnog_annot;
    io::CSVReader<EGGNOG_COL_NUM, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(path);
    while (in.read_row(qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                       best_og, cog_cat, eggnog_annot)) {
        // Check if the query matches one of our original transcriptome sequences
        QUERY_MAP_T::iterator it = (*mpQueryData->get_sequences_ptr()).find(qseqid);
        if (it != (*mpQueryData->get_sequences_ptr()).end()) {
            // EggNOG hit matches one of our original queries (from transcriptome)

            count_TOTAL_hits++;     // Increment number of EggNOG hits we got

            // Compile EggNOG results to be added to overall info
            EggnogResults = {};
            EggnogResults.seed_ortholog = seed_ortho;
            EggnogResults.seed_evalue = seed_e;
            EggnogResults.seed_score = seed_score;
            EggnogResults.predicted_gene = predicted_gene;
            EggnogResults.ogs = ogs;
            EggnogResults.raw_go = pFileSystem->list_to_vect(',', go_terms);    // Turn list into vect
            EggnogResults.raw_kegg = pFileSystem->list_to_vect(',', kegg);      // Turn list into vect
            get_tax_scope(tax_scope, EggnogResults);        // Map virNOG[6] to viridiplantae
            get_og_query(EggnogResults);               // Requires tax_scope first, pulls key to SQL database
            get_sql_data(EggnogResults, EGGNOG_DATABASE);
            EggnogResults.parsed_go = parse_go_list(go_terms,pEntapDatabase,',');

            it->second->set_eggnog_results(EggnogResults);  // Set EggNOG results to maintained data

            //  Analyze Gene Ontology Stats
            if (!EggnogResults.parsed_go.empty()) {
                count_total_go_hits++;
                it->second->QUERY_FLAG_SET(QuerySequence::QUERY_ONE_GO);
                for (auto &pair : EggnogResults.parsed_go) {
                    // pair - first: GO category, second; vector of terms
                    for (std::string &term : pair.second) {
                        // Count the terms we've found for individual category
                        go_combined_map[pair.first].add_value(term);
                        // Count the terms we've found overall (not category specific)
                        go_combined_map[GO_OVERALL_FLAG].add_value(term);
                    }
                }
            } else {
                count_no_go++;  // Increment number of EggNOG hits that did not have a GO term
            }

            // Analyze KEGG stats
            if (!kegg.empty()) {
                count_total_kegg_hits++;
                ct = (uint32) std::count(kegg.begin(), kegg.end(), ',');
                count_total_kegg_terms += ct + 1;
                it->second->QUERY_FLAG_SET(QuerySequence::QUERY_ONE_KEGG);
            } else {
                count_no_kegg++;
            }

            // Compile Taxonomic Orthogroup stats
            if (!EggnogResults.tax_scope_readable.empty()) {
                // Count number of unique taxonomic groups
                tax_scope_counter.add_value(EggnogResults.tax_scope_readable);
            }
        } else {
            // EggNOG hit does NOT match one of our original transcripts (must be some formatting error)
            throw ExceptionHandler("Sequence ID does not match in transcriptome: " + qseqid, ERR_ENTAP_PARSE_EGGNOG);
        }
    }
    // delete temp file
    pFileSystem->delete_file(path);

    // Close EggNOG SQL database
    EGGNOG_DATABASE.close();

    // Prepare output files
    out_no_hits_nucl = PATHS(mProcDir, OUT_UNANNOTATED_NUCL);
    out_no_hits_prot = PATHS(mProcDir, OUT_UNANNOTATED_PROT);
    out_hit_nucl     = PATHS(mProcDir, OUT_ANNOTATED_NUCL);
    out_hit_prot     = PATHS(mProcDir, OUT_ANNOTATED_PROT);
    std::ofstream file_no_hits_nucl(out_no_hits_nucl, std::ios::out | std::ios::app);
    std::ofstream file_no_hits_prot(out_no_hits_prot, std::ios::out | std::ios::app);
    std::ofstream file_hits_nucl(out_hit_nucl, std::ios::out | std::ios::app);
    std::ofstream file_hits_prot(out_hit_prot, std::ios::out | std::ios::app);

    FS_dprint("Success! Computing overall statistics...");
    // Find how many original sequences did/did not hit the EggNOG database
    for (auto &pair : *mpQueryData->get_sequences_ptr()) {
        if (!pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_EGGNOG_HIT)) {
            // Unannotated sequence
            if (!pair.second->get_sequence_n().empty()) file_no_hits_nucl << pair.second->get_sequence_n() << std::endl;
            if (!pair.second->get_sequence_p().empty()) file_no_hits_prot << pair.second->get_sequence_p() << std::endl;
            count_no_hits++;
        } else {
            // Annotated sequence
            if (!pair.second->get_sequence_n().empty()) file_hits_nucl << pair.second->get_sequence_n() << std::endl;
            if (!pair.second->get_sequence_p().empty()) file_hits_prot << pair.second->get_sequence_p() << std::endl;
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

    // Make sure we have hits before doing anything
    if (count_total_go_hits > 0) {

        //--------------------- Top Ten Taxonomic Scopes --------------//
        if (!tax_scope_counter.empty()) {
            // Setup graphing files
            std::string fig_txt_tax_bar = PATHS(mFigureDir, GRAPH_EGG_TAX_BAR_TXT);
            std::string fig_png_tax_bar = PATHS(mFigureDir, GRAPH_EGG_TAX_BAR_PNG);
            std::ofstream file_tax_bar(fig_txt_tax_bar, std::ios::out | std::ios::app);
            file_tax_bar << "Taxonomic Scope\tCount" << std::endl;

            ss << "\nTop 10 Taxonomic Scopes Assigned:";
            ct = 1;
            // Sort taxonomy scope
            tax_scope_counter.sort(true);
            for (auto &pair : tax_scope_counter._sorted) {
                if (ct > COUNT_TOP_TAX_SCOPE) break;
                percent = ((fp32)pair.second / tax_scope_counter._ct_total) * 100;
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
            mpGraphingManager->graph(graphingStruct);
        }
        //-------------------------------------------------------------//

        ss<<
          "\nTotal unique sequences with at least one GO term: " << count_total_go_hits <<
          "\nTotal unique sequences without GO terms: " << count_no_go <<
          "\nTotal GO terms assigned: " << go_combined_map[GO_OVERALL_FLAG]._ct_total;;

        for (uint16 lvl : mGoLevels) {
            for (auto &pair : go_combined_map) {
                if (pair.first.empty()) continue;
                // Count maps (biological/molecular/cellular/overall)
                fig_txt_go_bar = PATHS(mFigureDir, pair.first) + std::to_string(lvl)+GRAPH_GO_END_TXT;
                fig_png_go_bar = PATHS(mFigureDir, pair.first) + std::to_string(lvl)+GRAPH_GO_END_PNG;
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
                ss << "\nTotal "        << pair.first <<" terms (lvl="          << lvl << "): " << lvl_ct;
                ss << "\nTotal unique " << pair.first <<" terms (lvl="          << lvl << "): " << ct;
                ss << "\nTop 10 "       << pair.first <<" terms assigned (lvl=" << lvl << "): ";

                ct = 1;
                for (auto &pair2 : pair.second._sorted) {
                    if (ct > COUNT_TOP_GO) break;
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
                graphingStruct.software_flag = GRAPH_ONTOLOGY_FLAG;
                graphingStruct.graph_type = GRAPH_TOP_BAR_FLAG;
                mpGraphingManager->graph(graphingStruct);
            }
        }
        ss<<
          "\nTotal unique sequences with at least one pathway (KEGG) assignment: " << count_total_kegg_hits<<
          "\nTotal unique sequences without pathways (KEGG): " << count_no_kegg<<
          "\nTotal pathways (KEGG) assigned: " << count_total_kegg_terms;
    } else {
        // NO alignments for EggNOG, warn user!!. Don't fatal error
        FS_dprint("Warning: No EggNOG alignments found!!");
        ss << "\nWarning: No EggNOG Alignments";
    }
    out_msg = ss.str();
    pFileSystem->print_stats(out_msg);
    FS_dprint("Success! EggNOG results parsed");
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
    std::string                        blast;

    if (pFileSystem->file_exists(mOutHIts)) return;   // Should never occur here
    mBlastp ? blast = "blastp" : blast = "blastx";

    std::string                        annotation_base_flag;
    std::string                        annotation_std;
    std::string                        eggnog_command;
    std::string                        no_hit_out;

    annotation_base_flag = PATHS(_egg_out_dir, EGG_ANNOT_RESULTS);
    annotation_std       = PATHS(_egg_out_dir, EGG_ANNOT_STD);
    eggnog_command       = "python " + EGG_EMAPPER_EXE + " ";

    std::unordered_map<std::string,std::string> eggnog_command_map = {
            {"-i",mInFastaPath},
            {"--output",annotation_base_flag},
            {"--cpu",std::to_string(mThreads)},
            {"-m", "diamond"}
    };
    if (!mBlastp) eggnog_command_map["--translate"] = " ";
    if (pFileSystem->file_exists(mInFastaPath)) {
        for (auto &pair : eggnog_command_map)eggnog_command += pair.first + " " + pair.second + " ";
        if (TC_execute_cmd(eggnog_command, annotation_std) !=0) {
            pFileSystem->delete_file(annotation_base_flag);
            throw ExceptionHandler("Error executing eggnog mapper", ERR_ENTAP_RUN_ANNOTATION);
        }
        FS_dprint("Success! Results written to: " + annotation_base_flag);
    } else {
        throw ExceptionHandler("No input file found at: " + mInFastaPath, ERR_ENTAP_RUN_EGGNOG);
    }
    FS_dprint("Success! EggNOG execution complete");
}


// TODO remove
std::string ModEggnog::eggnog_format(std::string file) {
    FS_dprint("Reformatting EggNOG output: " + file);

    std::string out_path;
    std::string line;

    out_path = file + "_alt";
    pFileSystem->delete_file(out_path);
    std::ifstream in_file(file);
    std::ofstream out_file(out_path);
    while (getline(in_file,line)) {
        if (line.at(0) == '#' || line.empty()) continue;
        out_file << line << std::endl;
    }
    in_file.close();
    out_file.close();
    FS_dprint("Success! File printed to: " + out_path);
    return out_path;
}



/**
 * ======================================================================
 * Function bool ModEggnog::is_executable()
 *
 * Description          - Test command to see if eggnog is properly installed
 *                      - Called from UserInput.c
 *
 * Notes                - None
 *
 *
 * @return              - True/false successful execution or not
 * ======================================================================
 */
bool ModEggnog::is_executable() {
    std::string test_command;

    test_command = "python " +
            EGG_EMAPPER_EXE  +
            " --version";
    return TC_execute_cmd(test_command) == 0;
}

ModEggnog::~ModEggnog() {
    FS_dprint("Killing object - ModEggnog");
}


/**
 * ======================================================================
 * Function void ModEggnog::get_tax_scope(std::string &raw_scope,
                                    QuerySequence::EggnogResults &eggnogResults)
 *
 * Description          - Pulls the readable taxonomic scope from the EggnogLevels.h
 *                        file.
 *
 * Notes                - These may change, they are not pulled at configu (hardcoded)
 *
 * @param raw_scope     - Raw tax scope from EggNOG output (ie. virNOG[6])
 * @param eggnogResults - Current query sequence Eggnog struc
 *
 * @return              - None
 * ======================================================================
 */
void ModEggnog::get_tax_scope(std::string &raw_scope,
                                    QuerySequence::EggnogResults &eggnogResults) {
    // Lookup/Assign Tax Scope

    if (!raw_scope.empty()) {
        uint16 p = (uint16) (raw_scope.find("NOG"));
        if (p != std::string::npos) {
            eggnogResults.tax_scope = raw_scope.substr(0,p+3);
            eggnogResults.tax_scope_readable = EGGNOG_LEVELS[eggnogResults.tax_scope];
            return;
        }
    }
    eggnogResults.tax_scope  = raw_scope;
    eggnogResults.tax_scope_readable = "";
}


/**
 * ======================================================================
 * Function void ModEggnog::get_sql_data(QuerySequence::EggnogResults &eggnogResults,
 *                                          DatabaseHelper &database)
 *
 * Description          - Query EggNOG SQL database and pull relevant info
 *                      - Sets values in EggnogResults struct
 *
 * Notes                - None
 *
 * @param database      - DatabaseHelper object of Eggnog database
 * @param eggnogResults - Current query sequence Eggnog struc
 *
 * @return              - None
 * ======================================================================
 */
void ModEggnog::get_sql_data(QuerySequence::EggnogResults &eggnogResults, SQLDatabaseHelper &database) {
    // Lookup description, KEGG, protein domain from SQL database
    if (!eggnogResults.og_key.empty()) {
        std::vector<std::vector<std::string>>results;
        std::string sql_kegg;
        std::string sql_desc;
        std::string sql_protein;

        char *query = sqlite3_mprintf(
                "SELECT description, KEGG_freq, SMART_freq FROM og WHERE og=%Q",
                eggnogResults.og_key.c_str());
        try {
            results = database.query(query);
            sql_desc = results[0][0];
            sql_kegg = results[0][1];
            sql_protein = results[0][2];
            if (!sql_desc.empty() && sql_desc.find("[]") != 0) eggnogResults.description = sql_desc;
            if (!sql_kegg.empty() && sql_kegg.find("[]") != 0) {
                eggnogResults.sql_kegg = format_sql_data(sql_kegg);
            }
            if (!sql_protein.empty() && sql_protein.find("{}") != 0){
                eggnogResults.protein_domains = format_sql_data(sql_protein);
            }
        } catch (std::exception &e) {
            // Do not fatal error
            FS_dprint(e.what());
        }
    }
}


/**
 * ======================================================================
 * Function void ModEggnog::get_og_query(QuerySequence::EggnogResults &eggnogResults)
 *
 * Description          - Find specific OG that aligned to database (based on tax scope)
 *                      - Sets values in EggnogResults struct
 *
 * Notes                - None
 *
 * @param eggnogResults - Current query sequence Eggnog struc
 *
 * @return              - None
 * ======================================================================
 */
void ModEggnog::get_og_query(QuerySequence::EggnogResults &eggnogResults) {
    // Find OG query was assigned to
    std::string temp;
    if (!eggnogResults.ogs.empty()) {
        std::istringstream ss(eggnogResults.ogs);
        std::unordered_map<std::string,std::string> og_map; // Not fully used right now
        while (std::getline(ss,temp,',')) {
            uint16 p = (uint16) temp.find("@");
            og_map[temp.substr(p+1)] = temp.substr(0,p);
        }
        eggnogResults.og_key = "";
        if (og_map.find(eggnogResults.tax_scope) != og_map.end()) {
            eggnogResults.og_key = og_map[eggnogResults.tax_scope];
        }
    }
}


/**
 * ======================================================================
 * Function std::string ModEggnog::format_sql_data(std::string &input)
 *
 * Description          - Parses SQL data and puts it in generic format for EnTAP
 *
 * Notes                - None
 *
 * @param input         - Unformatted SQL data
 *
 * @return              - Formatted string
 * ======================================================================
 */
std::string ModEggnog::format_sql_data(std::string &input) {
    enum FLAGS {
        DOM    = 0x01,
        INNER  = 0x02,
        INNER2 = 0x04,
        STR    = 0x08,
        FOUND  = 0x10
    };

    unsigned char bracketFlag = 0x0;
    std::string output = "";

    for (char c : input) {
        if (c == '{') {
            bracketFlag |= DOM;
        } else if (c == '[' && !(bracketFlag & INNER)) {
            bracketFlag |= INNER;
            if (bracketFlag & DOM) output += " (";
        } else if (c == '[' && !(bracketFlag & INNER2)) {
            bracketFlag |= INNER2;
        } else if (c == '}') {
            bracketFlag &= ~DOM;
            bracketFlag &= ~FOUND;
            output = output.substr(0,output.length()-2); // Remove trailing ', '
        } else if (c == ']' && (bracketFlag & INNER2)) {
            bracketFlag &= ~INNER2;
            bracketFlag &= ~FOUND;
            output += ", ";
        } else if (c == ']' && (bracketFlag & INNER)) {
            bracketFlag &= ~INNER;
            bracketFlag &= ~FOUND;
            output = output.substr(0,output.length()-2); // Remove trailing ', '
            if (bracketFlag & DOM) output += "), ";
        } else if (c == '\"') {
            bracketFlag ^= STR;
            if (!(bracketFlag & STR)) bracketFlag |= FOUND;
            if (!(bracketFlag & INNER)) bracketFlag &= ~FOUND;
        } else {
            if (bracketFlag & FOUND) continue;
            if (bracketFlag & STR) output += c;
        }
    }
    return output;
}

/**
 * ======================================================================
 * Function bool ModEggnog::valid_input(boostPO::variables_map &)
 *
 * Description          - Ensure the input from the user is valid
 *                      - Called in UserInput
 *
 * Notes                - None
 *
 * @param user_map      - Map of user input (Boost)
 *
 * @return              - None
 * ======================================================================
 */
bool ModEggnog::valid_input(boostPO::variables_map &user_map) {
    return true;
}

ModEggnog::ModEggnog(std::string &out, std::string &in, std::string &ont,
                     bool blastp, std::vector<uint16> &lvls,EntapDataPtrs &entap_data,
                     std::string& eggnog_sql_path)
    :AbstractOntology(out, in,ont, blastp, lvls, entap_data){
    FS_dprint("Spawn Object - EggNOG");

    mEggnogDbPath = eggnog_sql_path;

    _egg_out_dir= PATHS(mOntologyDir, EGGNOG_DIRECTORY);
    mFigureDir = PATHS(_egg_out_dir, FIGURE_DIR);
    mProcDir   = PATHS(_egg_out_dir, PROCESSED_OUT_DIR);

    pFileSystem->delete_dir(mFigureDir);
    pFileSystem->delete_dir(mProcDir);

    pFileSystem->create_dir(_egg_out_dir);
    pFileSystem->create_dir(mFigureDir);
    pFileSystem->create_dir(mProcDir);
}

#endif
