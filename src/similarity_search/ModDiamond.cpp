/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
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
#include "ModDiamond.h"
#include "../QuerySequence.h"
#ifdef USE_BOOST
#include <boost/regex.hpp>
#else   // C++ libs
#include <regex>
#endif

std::vector<ENTAP_HEADERS> ModDiamond::DEFAULT_HEADERS = {
        ENTAP_HEADER_QUERY,
        ENTAP_HEADER_FRAME,
        ENTAP_HEADER_EXP_FPKM,
        ENTAP_HEADER_SIM_SUBJECT,
        ENTAP_HEADER_SIM_PERCENT,
        ENTAP_HEADER_SIM_ALIGN_LEN,
        ENTAP_HEADER_SIM_MISMATCH,
        ENTAP_HEADER_SIM_GAP_OPEN,
        ENTAP_HEADER_SIM_QUERY_S,
        ENTAP_HEADER_SIM_QUERY_E,
        ENTAP_HEADER_SIM_SUBJ_S,
        ENTAP_HEADER_SIM_SUBJ_E,
        ENTAP_HEADER_SIM_E_VAL,
        ENTAP_HEADER_SIM_COVERAGE,
        ENTAP_HEADER_SIM_TITLE,
        ENTAP_HEADER_SIM_SPECIES,
        ENTAP_HEADER_SIM_DATABASE,
        ENTAP_HEADER_SIM_CONTAM,
        ENTAP_HEADER_SIM_INFORM,
        ENTAP_HEADER_SIM_UNI_GO_BIO,
        ENTAP_HEADER_SIM_UNI_GO_CELL,
        ENTAP_HEADER_SIM_UNI_GO_MOLE,
        ENTAP_HEADER_SIM_UNI_KEGG,
        ENTAP_HEADER_SIM_UNI_DATA_XREF,
        ENTAP_HEADER_SIM_UNI_COMMENTS
};

ModDiamond::ModDiamond(std::string &execution_stage_path, std::string &in_hits, EntapDataPtrs &entap_data,
                       std::string &exe, vect_str_t &databases)
: AbstractSimilaritySearch(execution_stage_path, in_hits, entap_data, "DIAMOND", exe, databases){
    FS_dprint("Spawn Object - ModDiamond");

    _software_flag = SIM_DIAMOND;
}

EntapModule::ModVerifyData ModDiamond::verify_files() {
    // Transcriptome + database paths already verified
    ModVerifyData verify_data;
    std::string   database_name;        // Shortened name to be used for file naming
    std::string   file_name;            // File name that will be used for each database
    std::string   out_path;             // Full output path for each database alignment
    std::string   std_out;
    uint16 file_status = 0;

    verify_data.files_exist = true;

    for (std::string &data_path : _database_paths) {
        FS_dprint("Verifying previous execution of database: " + data_path + "...");

        database_name = get_database_shortname(data_path);

        // set output file name (hits will be printed to this for diamond)
        file_name = get_output_path(data_path);

        // set full path output for this database
        out_path = PATHS(_mod_out_dir, file_name);
        std_out  = PATHS(_mod_out_dir, file_name) + FileSystem::EXT_STD;

        // add mapping of output file to shortened database name
        _path_to_database[out_path] = database_name;

        // Check if file exists/can be read/empty
        file_status = _pFileSystem->get_file_status(out_path);
        if (file_status != 0) {
            FS_dprint("File for database " + database_name + " does not exist.\n" + out_path);
            // If we need to execute against ANY database
            verify_data.files_exist = false;
            // delete file just in case it is corrupt/empty
            _pFileSystem->delete_file(out_path);
        } else {
            // File found + is 'legit', can skip execution for it
            FS_dprint("File for database " + database_name + " exists, skipping...\n" + out_path);
        }

        verify_data.output_paths.push_back(out_path);
        _output_paths.push_back(out_path);
    }

    FS_dprint("Success! Verified files for DIAMOND, continuing...");

    return verify_data;
}

bool ModDiamond::is_executable(std::string &exe) {
    std::string test_command;
    TerminalData terminalData;

    test_command = exe + " --version";

    terminalData.command = test_command;
    terminalData.print_files = false;

    return TC_execute_cmd(terminalData) == 0;
}

void ModDiamond::execute() {
    std::string output_path;
    uint16 file_status = 0;
    SimSearchCmd simSearchCmd;

    FS_dprint("Executing DIAMOND for necessary files....");

    for (std::string &database_path : _database_paths) {
        output_path = get_output_path(database_path);

        file_status = _pFileSystem->get_file_status(output_path);
        if (file_status != 0) {
            // If file does not exist or cannot be read, execute diamond
            FS_dprint("File not found, executing against database at: " + database_path);

            simSearchCmd = {};
            simSearchCmd.database_path = database_path;
            simSearchCmd.output_path   = output_path;
            simSearchCmd.std_out_path  = output_path + FileSystem::EXT_STD;
            simSearchCmd.threads       = (uint16)_threads;
            simSearchCmd.query_path    = _in_hits;
            simSearchCmd.eval          = _e_val;
            simSearchCmd.tcoverage     = _tcoverage;
            simSearchCmd.qcoverage     = _qcoverage;
            simSearchCmd.exe_path      = _exe_path;

            try {
                run_blast(&simSearchCmd, true);
            } catch (const ExceptionHandler &e ){
                throw e;
            }

            FS_dprint("Success! Results written to: " + output_path);
        }
    }
}

bool ModDiamond::run_blast(AbstractSimilaritySearch::SimSearchCmd *cmd, bool use_defaults) {
    std::string     diamond_cmd;
    TerminalData    terminalData;
    int32           err_code;
    bool            ret = true;

    diamond_cmd = cmd->exe_path + " ";

    if (cmd->blastp) {
        diamond_cmd += BLASTP_STR;
    } else {
        diamond_cmd += BLASTX_STR;
    }

    diamond_cmd += " -d " + cmd->database_path;
    diamond_cmd += " --query-cover " + std::to_string(cmd->qcoverage);
    diamond_cmd += " --subject-cover " + std::to_string(cmd->tcoverage);
    diamond_cmd += " --evalue " + std::to_string(cmd->eval);

    diamond_cmd += " --more-sensitive --top 3";

    diamond_cmd += " -q " + cmd->query_path;
    diamond_cmd += " -o " + cmd->output_path;
    diamond_cmd += " -p " + std::to_string(cmd->threads);

    terminalData.command        = diamond_cmd;
    terminalData.base_std_path  = cmd->std_out_path;
    terminalData.print_files    = true;

    err_code = TC_execute_cmd(terminalData);

    // will change at some point
    if (err_code != 0) {
        // delete output file if run failed
        _pFileSystem->delete_file(cmd->output_path);
        throw ExceptionHandler("Error with database located at: " + cmd->database_path + "\nDIAMOND Error: " +
            terminalData.err_stream, ERR_ENTAP_RUN_SIM_SEARCH_RUN);
    }

    return ret;
}

void ModDiamond::parse() {
    bool                is_uniprot;
    uint32              uniprot_attempts=0;
    uint16              file_status=0;
    std::string         database_shortname;
    std::string         species;
    QuerySequence::SimSearchResults simSearchResults;
    TaxEntry            taxEntry;
    std::pair<bool, std::string> contam_info;

    // ------------------ Read from DIAMOND output ---------------------- //
    std::string qseqid, sseqid, stitle, database_name,pident, bitscore,
            length, mismatch, gapopen, qstart, qend, sstart, send;
    fp64  evalue, coverage;

    // ----------------------------------------------------------------- //

    FS_dprint("Beginning to filter individual DIAMOND files...");

    // disable UniProt headers until we know we have a hit
    EM_set_uniprot_headers(false);

    for (std::string &output_path : _output_paths) {
        FS_dprint("DIAMOND file located at " + output_path + " being parsed");

        // reset uniprot info for each database
        is_uniprot = false;
        uniprot_attempts = 0;
        simSearchResults = {};

        // ensure file exists
        file_status = _pFileSystem->get_file_status(output_path);
        if (file_status != 0) {
            throw ExceptionHandler("File not found or empty: " + output_path, ERR_ENTAP_RUN_SIM_SEARCH_FILTER);
        }

        // setup individual database directories for stats/figures
        database_shortname = _path_to_database[output_path];

        // Begin using CSVReader lib to parse data
        io::CSVReader<DMND_COL_NUMBER, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(output_path);
        while (in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                           qstart, qend, sstart, send, evalue, bitscore, coverage,stitle)) {
            simSearchResults = {};

            // Get pointer to sequence in overall map
            QuerySequence *query = _pQUERY_DATA->get_sequence(qseqid);
            if (query == nullptr) {
                throw ExceptionHandler("Unable to find sequence in transcriptome: " + qseqid + " from file: " + output_path,
                                       ERR_ENTAP_RUN_SIM_SEARCH_FILTER);
            }

            // get species from database alignment (using boost regex for now)
            species = get_species(stitle);
            // get taxonomic information with species
            taxEntry = _pEntapDatabase->get_tax_entry(species);
            // get contaminant information
            contam_info = is_contaminant(taxEntry.lineage, _contaminants);

            // Check if this is a UniProt match and pull back info if so
            if (is_uniprot) {
                // Get uniprot info
                is_uniprot_entry(sseqid, simSearchResults.uniprot_info);
            } else {
                if (uniprot_attempts <= UNIPROT_ATTEMPTS) {
                    // First UniProt match assumes the rest are UniProt as well in database
                    is_uniprot = is_uniprot_entry(sseqid, simSearchResults.uniprot_info);
                    if (!is_uniprot) {
                        uniprot_attempts++;
                    } else {
                        FS_dprint("Database file at " + output_path + "\nDetermined to be UniProt");
                        _pQUERY_DATA->DATA_FLAG_SET(QueryData::UNIPROT_MATCH);
                        EM_set_uniprot_headers(true);
                    }
                } // Else, database is NOT UniProt after # of attempts
            }

            // Compile sim search data
            simSearchResults.database_path = output_path;
            simSearchResults.qseqid = qseqid;
            simSearchResults.sseqid = sseqid;
            simSearchResults.pident = pident;
            simSearchResults.length = length;
            simSearchResults.mismatch = mismatch;
            simSearchResults.gapopen = gapopen;
            simSearchResults.qstart = qstart;
            simSearchResults.qend = qend;
            simSearchResults.sstart = sstart;
            simSearchResults.send = send;
            simSearchResults.stitle = stitle;
            simSearchResults.bit_score = bitscore;
            simSearchResults.lineage = taxEntry.lineage;
            simSearchResults.species = species;
            simSearchResults.e_val_raw = evalue;
            simSearchResults.e_val = float_to_sci(evalue,2);
            simSearchResults.coverage_raw = coverage;
            simSearchResults.coverage = float_to_string(coverage);
            simSearchResults.contaminant = contam_info.first;
            simSearchResults.contam_type = contam_info.second;
            simSearchResults.contaminant ? simSearchResults.yes_no_contam = YES_FLAG :
                    simSearchResults.yes_no_contam  = NO_FLAG;
            simSearchResults.is_informative = is_informative(stitle, _uninformative_vect);
            simSearchResults.is_informative ? simSearchResults.yes_no_inform = YES_FLAG :
                    simSearchResults.yes_no_inform  = NO_FLAG;

            query->add_alignment(_execution_state, _software_flag,
                    simSearchResults, output_path, _input_lineage);
        } // END WHILE LOOP

        // Finished parsing and adding to alignment data, being to calc stats
        FS_dprint("File parsed, calculating statistics and writing output...");
        calculate_best_stats(false,output_path);
        FS_dprint("Success!");
    } // END FOR LOOP

    FS_dprint("Calculating overall Similarity Searching statistics...");
    calculate_best_stats(true);
    FS_dprint("Success!");
}

typedef std::map<std::string,std::map<std::string,uint32>> graph_sum_t;

void ModDiamond::calculate_best_stats (bool is_final, std::string database_path) {

    GraphingData                graphingStruct;
    std::string                 species;
    std::string                 database_shortname;
    std::string                 figure_base;
    std::string                 frame;
    std::string                 base_path;
    std::string                 temp_file_path;
    std::string                 contam;
    std::stringstream           ss;
    uint64                      count_no_hit=0;
    uint64                      count_contam=0;
    uint64                      count_filtered=0;
    uint64                      count_informative=0;
    uint64                      count_uninformative=0;
    uint64                      count_unselected=0;
    uint64                      count_TOTAL_alignments=0;
    uint32                      ct;
    fp64                        percent;
    fp64                        contam_percent;
    Compair<std::string>        contam_counter;
    Compair<std::string>        species_counter;
    Compair<std::string>        contam_species_counter;
    graph_sum_t                 graphing_sum_map;

    // Set up output directories (processed directory cleared earlier so these will be empty)
    if (is_final) {
        // Overall results across databases
        base_path = _overall_results_dir;
        database_shortname = "";
    } else {
        // Individual database results
        database_shortname = _path_to_database[database_path];
        base_path   = PATHS(_proc_dir, database_shortname);
    }
    figure_base = PATHS(base_path, FIGURE_DIR);
    _pFileSystem->create_dir(base_path);
    _pFileSystem->create_dir(figure_base);

    // Open contam best hit tsv file and print headers
    std::string out_best_contams_filepath = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_HITS_CONTAM);
    _pQUERY_DATA->start_alignment_files(out_best_contams_filepath, DEFAULT_HEADERS, 0, nullptr);

    // Open best hits files
    std::string out_best_hits_filepath = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_HITS);
    _pQUERY_DATA->start_alignment_files(out_best_hits_filepath, DEFAULT_HEADERS, 0, nullptr);

    // Open best hits files with no contaminants
    std::string out_best_hits_no_contams = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_HITS_NO_CONTAM);
    _pQUERY_DATA->start_alignment_files(out_best_hits_no_contams, DEFAULT_HEADERS, 0, nullptr);

    // Open unselected hits, so every hit that was not the best hit (tsv)
    std::string out_unselected_tsv  = PATHS(base_path, SIM_SEARCH_DATABASE_UNSELECTED + FileSystem::EXT_TSV);
    std::ofstream file_unselected_hits(out_unselected_tsv, std::ios::out | std::ios::app);

    // Open no hits file (fasta nucleotide)
    std::string out_no_hits_fa_nucl = PATHS(base_path, SIM_SEARCH_DATABASE_NO_HITS + FileSystem::EXT_FNN);
    std::ofstream file_no_hits_nucl(out_no_hits_fa_nucl, std::ios::out | std::ios::app);

    // Open no hits file (fasta protein)
    std::string out_no_hits_fa_prot  = PATHS(base_path, SIM_SEARCH_DATABASE_NO_HITS + FileSystem::EXT_FAA);
    std::ofstream file_no_hits_prot(out_no_hits_fa_prot, std::ios::out | std::ios::app);

    // ------------------- Setup graphing files ------------------------- //

    std::string graph_species_txt_path           = PATHS(figure_base, GRAPH_SPECIES_BAR_TXT);
    std::string graph_species_png_path           = PATHS(figure_base, GRAPH_SPECIES_BAR_PNG);
    std::string graph_contam_txt_path            = PATHS(figure_base, GRAPH_CONTAM_BAR_TXT);
    std::string graph_contam_png_path            = PATHS(figure_base, GRAPH_CONTAM_BAR_PNG);
    std::string graph_sum_txt_path               = PATHS(figure_base, GRAPH_DATABASE_SUM_TXT);
    std::string graph_sum_png_path               = PATHS(figure_base, GRAPH_DATABASE_SUM_PNG);

    std::ofstream graph_species_file(graph_species_txt_path, std::ios::out | std::ios::app);
    std::ofstream graph_contam_file(graph_contam_txt_path, std::ios::out | std::ios::app);
    std::ofstream graph_sum_file(graph_sum_txt_path, std::ios::out | std::ios::app);

    // ------------------------------------------------------------------ //

    // Print headers to relevant tsv files
    _pFileSystem->print_headers(file_unselected_hits, DEFAULT_HEADERS, FileSystem::DELIM_TSV);


    try {
        graph_species_file << "Species\tCount"     << std::endl;
        graph_contam_file  << "Contaminant Species\tCount" << std::endl;
        graph_sum_file     << "Category\tCount"    << std::endl;

        // Cycle through all sequences
        for (auto &pair : *_pQUERY_DATA->get_sequences_ptr()) {
            // Check if original sequences have hit a database
            if (!pair.second->hit_database(SIMILARITY_SEARCH, SIM_DIAMOND, database_path)) {
                // Did NOT hit a database during sim search
                // Do NOT log if it was never blasted
                if ((pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN) && _blastp) ||
                    (!pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN) && !_blastp)) {
                    // Protein/nucleotide did not hit database
                    count_no_hit++;
                    file_no_hits_nucl << pair.second->get_sequence_n() << std::endl;
                    file_no_hits_prot << pair.second->get_sequence_p() << std::endl;
                    // Graphing
                    frame = pair.second->getFrame();
                    if (graphing_sum_map[frame].find(NO_HIT_FLAG) != graphing_sum_map[frame].end()) {
                        graphing_sum_map[frame][NO_HIT_FLAG]++;
                    } else graphing_sum_map[frame][NO_HIT_FLAG] = 1;
                } else {
                    pair.second->QUERY_FLAG_SET(QuerySequence::QUERY_BLASTED);
                }
            } else {
                // HIT a database during sim search

                QuerySequence::SimSearchResults *sim_search_data;
                QuerySequence::SimSearchAlignment *best_hit;
                // Process unselected hits for non-final analysis and set best hit pointer
                if (is_final) {
                    best_hit =
                            pair.second->get_best_hit_alignment<QuerySequence::SimSearchAlignment>(
                                    SIMILARITY_SEARCH, SIM_DIAMOND,"");
                    sim_search_data = best_hit->get_results();
                } else {
                    best_hit = pair.second->get_best_hit_alignment<QuerySequence::SimSearchAlignment>(
                            SIMILARITY_SEARCH, SIM_DIAMOND,database_path);
                    QuerySequence::align_database_hits_t *alignment_data =
                            pair.second->get_database_hits(database_path,SIMILARITY_SEARCH, SIM_DIAMOND);
                    sim_search_data = best_hit->get_results();
                    for (auto &hit : *alignment_data) {
                        count_TOTAL_alignments++;
                        if (hit != best_hit) {  // If this hit is not the best hit
                            file_unselected_hits << hit->print_delim(DEFAULT_HEADERS, 0, FileSystem::DELIM_TSV) << std::endl;
                            count_unselected++;
                        } else {
                            ;   // Do notthing
                        }
                    }
                }
                count_filtered++;   // increment best hit

                // Write to best hits files
                _pQUERY_DATA->add_alignment_data(out_best_hits_filepath, DEFAULT_HEADERS, pair.second, 0);

                frame = pair.second->getFrame();     // Used for graphing
                species = sim_search_data->species;

                // Determine contaminant information and print to files
                if (sim_search_data->contaminant) {
                    // Species is considered a contaminant
                    count_contam++;
                    _pQUERY_DATA->add_alignment_data(out_best_contams_filepath, DEFAULT_HEADERS, pair.second, 0);

                    contam = sim_search_data->contam_type;
                    contam_counter.add_value(contam);
                    contam_species_counter.add_value(species);
                } else {
                    // Species is NOT a contaminant, print to files
                    _pQUERY_DATA->add_alignment_data(out_best_hits_no_contams, DEFAULT_HEADERS, pair.second, 0);
                }

                // Count species type
                species_counter.add_value(species);

                // Check if this is an informative alignment and respond accordingly
                if (sim_search_data->is_informative) {
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
            }
        }
    } catch (const std::exception &e){throw ExceptionHandler(e.what(), ERR_ENTAP_RUN_SIM_SEARCH_FILTER);}

    try {
        _pQUERY_DATA->end_alignment_files(out_best_contams_filepath);
        _pQUERY_DATA->end_alignment_files(out_best_hits_filepath);
        _pQUERY_DATA->end_alignment_files(out_best_hits_no_contams);

        _pFileSystem->close_file(file_no_hits_nucl);
        _pFileSystem->close_file(file_no_hits_prot);
        _pFileSystem->close_file(file_unselected_hits);
    } catch (const ExceptionHandler &e) {throw e;}

    // ------------ Calculate statistics and print to output ------------ //
    ss<<std::fixed<<std::setprecision(2);

    // Different headers if final analysis or database specific analysis
    if (is_final) {
        _pFileSystem->format_stat_stream(ss, "Compiled Similarity Search - DIAMOND - Best Overall");
    } else {
        _pFileSystem->format_stat_stream(ss, "Similarity Search - DIAMOND - " + database_shortname);
        ss <<
           "Search results:\n"            << database_path <<
           "\n\tTotal alignments: "               << count_TOTAL_alignments   <<
           "\n\tTotal unselected results: "       << count_unselected      <<
           "\n\t\tWritten to: "                   << out_unselected_tsv;
    }

    // If overall alignments are 0, then throw error
    if (is_final && count_filtered == 0) {
        throw ExceptionHandler("No alignments found during Similarity Searching!",
                               ERR_ENTAP_RUN_SIM_SEARCH_FILTER);
    }

    // If no total or filealignments for this database, return and warn user
    if (!is_final && (count_TOTAL_alignments == 0 || count_filtered == 0)) {
        ss << "WARNING: No alignments for this database";
        std::string out_msg = ss.str() + "\n";
        _pFileSystem->print_stats(out_msg);
        return;
    }

    // Sort counters
    contam_species_counter.sort(true);
    species_counter.sort(true);

    contam_percent = ((fp64) count_contam / count_filtered) * 100;

    ss <<
       "\n\tTotal unique transcripts with an alignment: " << count_filtered <<
       "\n\t\tReference transcriptome sequences with an alignment (FASTA):\n\t\t\t" << out_best_hits_filepath <<
       "\n\t\tSearch results (TSV):\n\t\t\t" << out_best_hits_filepath <<
       "\n\tTotal unique transcripts without an alignment: " << count_no_hit <<
       "\n\t\tReference transcriptome sequences without an alignment (FASTA):\n\t\t\t" << out_no_hits_fa_prot;
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
            graph_sum_file << pair.first << "\t" << INFORMATIVE_FLAG << "\t" << pair.second[INFORMATIVE_FLAG]
                           << "\n";
        }
    }
    ss <<
       "\n\tTotal unique uninformative alignments: " << count_uninformative;
    if (graphing_sum_map.size() > 1) {
        for (auto &pair : graphing_sum_map) {
            // Frame -> Map of uninform/inform/no hits
            ss << "\n\t\t" << pair.first << "(" << pair.second[UNINFORMATIVE_FLAG] << ")";
            graph_sum_file << pair.first << "\t" << UNINFORMATIVE_FLAG << "\t" << pair.second[UNINFORMATIVE_FLAG]
                           << "\n";
        }
    }

    ss <<
       "\n\tTotal unique contaminants: " << count_contam <<
       "(" << contam_percent << "%): " <<
       "\n\t\tTranscriptome reference sequences labeled as a contaminant (FASTA):\n\t\t\t"
       << out_best_contams_filepath <<
       "\n\t\tTranscriptome reference sequences labeled as a contaminant (TSV):\n\t\t\t" << out_best_contams_filepath;


    // ********** Contaminant Calculations ************** //
    if (count_contam > 0) {
        ss << "\n\t\tFlagged contaminants (all % based on total contaminants):";
        for (auto &pair : contam_counter._data) {
            percent = ((fp64) pair.second / count_contam) * 100;
            ss
                    << "\n\t\t\t" << pair.first << ": " << pair.second << "(" << percent << "%)";
        }
        ss << "\n\t\tTop " << COUNT_TOP_SPECIES << " contaminants by species:";
        ct = 1;
        for (auto &pair : contam_species_counter._sorted) {
            if (ct > COUNT_TOP_SPECIES) break;
            percent = ((fp64) pair.second / count_contam) * 100;
            ss
                    << "\n\t\t\t" << ct << ")" << pair.first << ": "
                    << pair.second << "(" << percent << "%)";
            graph_contam_file << pair.first << '\t' << std::to_string(pair.second) << std::endl;
            ct++;
        }
    }

    ss << "\n\tTop " << COUNT_TOP_SPECIES << " alignments by species:";
    ct = 1;
    for (auto &pair : species_counter._sorted) {
        if (ct > COUNT_TOP_SPECIES) break;
        percent = ((fp64) pair.second / count_filtered) * 100;
        ss
                << "\n\t\t\t" << ct << ")" << pair.first << ": "
                << pair.second << "(" << percent << "%)";
        graph_species_file << pair.first << '\t' << std::to_string(pair.second) << std::endl;
        ct++;
    }
    std::string out_msg = ss.str() + "\n";
    _pFileSystem->print_stats(out_msg);


    // ------------------------------------------------------------------ //
    // ********* Graphing Handle ********** //
    graphingStruct.software_flag = GRAPH_SOFTWARE_FLAG;
    _pFileSystem->close_file(graph_contam_file);
    _pFileSystem->close_file(graph_species_file);
    _pFileSystem->close_file(graph_sum_file);
    if (count_contam > 0) {
        graphingStruct.fig_out_path   = graph_contam_png_path;
        graphingStruct.graph_title    = database_shortname + GRAPH_CONTAM_TITLE;
        graphingStruct.text_file_path = graph_contam_txt_path;
        graphingStruct.graph_type     = GRAPH_BAR_FLAG;
        _pGraphingManager->graph(graphingStruct);
    }
    graphingStruct.fig_out_path   = graph_species_png_path;
    graphingStruct.graph_title    = database_shortname + GRAPH_SPECIES_TITLE;
    graphingStruct.text_file_path = graph_species_txt_path;
    graphingStruct.graph_type     = GRAPH_BAR_FLAG;
    _pGraphingManager->graph(graphingStruct);

    graphingStruct.fig_out_path   = graph_sum_png_path;
    graphingStruct.graph_title    = database_shortname + GRAPH_DATABASE_SUM_TITLE;
    graphingStruct.text_file_path = graph_sum_txt_path;
    graphingStruct.graph_type     = GRAPH_SUM_FLAG;
    _pGraphingManager->graph(graphingStruct);

    // check if final - different graph
    // ************************************ //
}
