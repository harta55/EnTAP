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
#include "ModDiamond.h"
#include "../QuerySequence.h"
#include "../QueryAlignment.h"
#include "../QueryData.h"
#include "../GraphingManager.h"
#include "../ExceptionHandler.h"

#ifdef USE_BOOST
#include <boost/regex.hpp>
#else   // C++ libs
#include <regex>
#endif

// Test defines
//#define TEST_SIM_PARSE_01

std::vector<ENTAP_HEADERS> ModDiamond::DEFAULT_HEADERS = {
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
        ENTAP_HEADER_SIM_TAXONOMIC_LINEAGE,
        ENTAP_HEADER_SIM_DATABASE,
        ENTAP_HEADER_SIM_CONTAM,
        ENTAP_HEADER_SIM_INFORM
};

std::vector<ENTAP_HEADERS> ModDiamond::UNIPROT_HEADERS = {
        ENTAP_HEADER_SIM_UNI_DATA_XREF,
        ENTAP_HEADER_SIM_UNI_COMMENTS,
        ENTAP_HEADER_SIM_UNI_KEGG,
        ENTAP_HEADER_SIM_UNI_GO_BIO,
        ENTAP_HEADER_SIM_UNI_GO_CELL,
        ENTAP_HEADER_SIM_UNI_GO_MOLE
};

/**
 * ======================================================================
 * Function void ModDiamond::ModDiamond(std::string &execution_stage_path,
 *                     std::string &fasta_path, EntapDataPtrs &entap_data,
                       std::string &exe, vect_str_t &databases)
 *
 * Description          - ModDiamond constructor, calls super to init
 *                        paths/set member variables
 *
 * Notes                - Constructor
 *
 * @param execution_stage_path - Absolute directory path to execution stage (Similarity Search)
 * @param in_hits              - Absolute path to FASTA transcriptome
 * @param entap_data           - Structure of necessary pointers for execution
 * @param exe                  - Module execution method (DIAMOND exe file)
 * @param databases            - Vector of databases to search against
 *
 *
 * @return              - ModDiamond object
 *
 * =====================================================================
 */
ModDiamond::ModDiamond(std::string &execution_stage_path, std::string &fasta_path, EntapDataPtrs &entap_data,
                       vect_str_t &databases)
: AbstractSimilaritySearch(execution_stage_path, fasta_path, entap_data, "DIAMOND", DEFAULT_HEADERS, databases){

    FS_dprint("Spawn Object - ModDiamond");
    mParsedFile = false;
    mSoftwareFlag = SIM_DIAMOND;
    mExePath = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_DIAMOND_EXE);
    mpFileSystem->delete_dir(mFigureDir);   // Don't need for DIAMOND, may change
}

/**
 * ======================================================================
 * Function EntapModule::ModVerifyData ModDiamond::verify_files()
 *
 * Description          - Determines whether DIAMOND has already been ran
 *                        with the same input transcriptome
 *
 * Notes                - None
 *
 *
 * @return              - Structure of ModVerifyData if files exist already
 *
 * =====================================================================
 */
EntapModule::ModVerifyData ModDiamond::verify_files() {
    // Transcriptome + database paths already verified
    ModVerifyData verify_data;
    std::string   database_name;        // Shortened name to be used for file naming
    std::string   out_path;             // Full output path for each database alignment
    uint16 file_status = 0;             // File statuses (empty, doesn't exist, etc...)

    verify_data.files_exist = true;

    for (std::string &data_path : mDatabasePaths) {
        FS_dprint("Verifying previous execution of database: " + data_path + "...");

        database_name = get_database_shortname(data_path);

        // set full path output for this database
        out_path = get_database_output_path(data_path);

        // add mapping of output file to shortened database name
        mPathToDatabase[out_path] = database_name;

        // Check if file exists/can be read/empty
        file_status = mpFileSystem->get_file_status(out_path);
        if (file_status != 0) {
            FS_dprint("File for database " + database_name + " does not exist.\n" + out_path);
            // If we need to execute against ANY database
            verify_data.files_exist = false;
            // delete file just in case it is corrupt/empty
            mpFileSystem->delete_file(out_path);
        } else {
            // File found + is 'legit', can skip execution for it
            FS_dprint("File for database " + database_name + " exists, skipping...\n" + out_path);
        }

        verify_data.output_paths.push_back(out_path);   // Add paths to verify data (not currently used)
        mOutputPaths.push_back(out_path);
    }

    FS_dprint("Success! Verified files for DIAMOND, continuing...");

    return verify_data;
}

/**
 * ======================================================================
 * Function bool ModDiamond::is_executable()
 *
 * Description          - Determines whether DIAMOND can execute on current system
 *
 * Notes                - None
 *
 *
 * @return              - BOOL is execution is possible
 *
 * =====================================================================
 */
bool ModDiamond::is_executable(std::string &exe) {
    TerminalData terminalData;  // Terminal data structure

    terminalData.command = exe + " --version";
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;

    return TC_execute_cmd(terminalData) == 0;
}

/**
 * ======================================================================
 * Function void ModDiamond::execute()
 *
 * Description          - Generates command for DIAMOND execution, calls
 *                        run_blast to execute DIAMOND
 *
 * Notes                - Databases must be configured for DIAMOND
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModDiamond::execute() {
    std::string output_path;    // Absolute path to output file from DIAMOND execution
    uint16 file_status = 0;     // File statuses (empty, cannot be read, etc...)
    SimSearchCmd simSearchCmd;  // DIAMOND commands

    FS_dprint("Executing DIAMOND for necessary files....");

    for (std::string &database_path : mDatabasePaths) {
        output_path = get_database_output_path(database_path);

        file_status = mpFileSystem->get_file_status(output_path);
        if (file_status != 0) {
            // If file does not exist or cannot be read, execute diamond
            FS_dprint("File not found, executing against database at: " + database_path);

            simSearchCmd = {};
            simSearchCmd.database_path = database_path;
            simSearchCmd.output_path   = output_path;
            simSearchCmd.std_out_path  = output_path + FileSystem::EXT_STD;
            simSearchCmd.threads       = (uint16)mThreads;
            simSearchCmd.query_path    = mInputTranscriptome;
            simSearchCmd.eval          = mEVal;
            simSearchCmd.tcoverage     = mTCoverage;
            simSearchCmd.qcoverage     = mQCoverage;
            simSearchCmd.exe_path      = mExePath;
            simSearchCmd.blastp        = mBlastp;

            try {
                run_blast(&simSearchCmd, true);
            } catch (const ExceptionHandler &e ){
                throw e;
            }

            FS_dprint("Success! Results written to: " + output_path);
        }
    }
}

/**
 * ======================================================================
 * Function bool ModDiamond::run_blast()
 *
 * Description          - Executes DIAMOND with provided command or uses defaults
 *                      - If defaults selected, only DIAMOND database/blast
 *                        command will be read
 *
 * Notes                - Databases must be configured for DIAMOND
 *
 * @param cmd           - DIAMOND flags for blast
 * @param use_defaults  - BOOL to use default values. If TRUE, override with defaults
 *
 * @return              - BOOL for successful execution
 *
 * =====================================================================
 */
bool ModDiamond::run_blast(AbstractSimilaritySearch::SimSearchCmd *cmd, bool use_defaults) {
    TerminalData    terminalData;   // Terminal data
    command_map_t   tc_commands;    // Terminal command map
    int32           err_code;       // Error codee from terminal execution
    bool            ret = true;     // Return value, if execution has succeeded
    std::string     temp_exe;       // DIAMOND needs blastp/x directly after DIAMOND exe

    // Overwrite values if we should use defaults
    if (use_defaults) {
        cmd->qcoverage = mQCoverage;
        cmd->tcoverage = mTCoverage;
        cmd->threads = (uint16)mThreads;
    }

    if (cmd->blastp) {
        temp_exe = cmd->exe_path + " " + CMD_BLASTP;
    } else {
        temp_exe = cmd->exe_path + " " + CMD_BLASTX;
    }

    tc_commands.emplace(CMD_DATABASE, cmd->database_path);
    tc_commands.emplace(CMD_QUERY_PATH, cmd->query_path);
    tc_commands.emplace(CMD_QUERY_COVERAGE, std::to_string(cmd->qcoverage));
    tc_commands.emplace(CMD_SUBJECT_COVERAGE, std::to_string(cmd->tcoverage));
    tc_commands.emplace(CMD_EVALUE, std::to_string(cmd->eval));
    tc_commands.emplace(CMD_MORE_SENSITIVE, TC_NULL_ARGUMENT);
    tc_commands.emplace(CMD_TOP_ALIGNMENTS, std::to_string(CMD_DEFAULT_TOP_ALIGN));

    tc_commands.emplace(CMD_OUTPUT_PATH, cmd->output_path);
    tc_commands.emplace(CMD_THREADS, std::to_string(cmd->threads));
    tc_commands.emplace(CMD_OUTPUT_FORMAT, CMD_DEFAULT_OUTPUT_FORMAT);

    terminalData.command        = TC_generate_command(tc_commands, temp_exe);
    terminalData.base_std_path  = cmd->std_out_path;
    terminalData.print_files    = true;
    terminalData.suppress_std_err = false;

    err_code = TC_execute_cmd(terminalData);

    // will change at some point
    if (err_code != 0) {
        // delete output file if run failed
        mpFileSystem->delete_file(cmd->output_path);
        throw ExceptionHandler("Error with database located at: " + cmd->database_path + "\nDIAMOND Error: " +
            terminalData.err_stream, ERR_ENTAP_RUN_SIM_SEARCH_RUN);
    }
    return ret;
}

/**
 * ======================================================================
 * Function void ModDiamond::parse()
 *
 * Description          - Parses and compiles data from DIAMOND output
 *                      - Adds this to QueryData class which determines best hits
 *
 * Notes                - None
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModDiamond::parse() {
    bool                is_uniprot;                     // BOOL is current database a UniProt database
    uint32              uniprot_attempts=0;             // Number of attempts to determine if database is UniProt
    uint16              file_status=0;                  // File statuses (defined FileSystem.h)
    std::string         species;                        // Species from subject sequence
    QuerySequence::SimSearchResults simSearchResults;   // Compiled similarity search results
    TaxEntry            taxEntry;                       // Entry from Tarxonomic database
    std::pair<bool, std::string> contam_info;           // Contaminate information
    std::stringstream   ss;                             // Output string stream

    // ------------------ Read from DIAMOND output ----------------------- //
    std::string qseqid;             // Sequence ID of query sequence
    std::string sseqid;             // Sequence ID of subject/target sequence (from database)
    std::string stitle;             // Title of subject sequence pulled from database
    std::string pident;             // Percent identical
    std::string bitscore;           // Bit score
    std::string length;             // Length (bp)
    std::string mismatch;           // Mismatches
    std::string gapopen;            // Gap open
    std::string qstart;             // Query start bp
    std::string qend;               // Query end bp
    std::string sstart;             // Subject start bp
    std::string send;               // SUbject end bp
    fp64 evalue;                    // E-value
    fp64 coverage;                  // Coverage
    // ------------------------------------------------------------------ //

#ifdef TEST_SIM_PARSE_01
    uint16 test_ct = 0;
    const uint16 TEST_COUNT_MAX = 100;       // Parse 100 lines only of each file
#endif

    FS_dprint("Beginning to filter individual DIAMOND files...");

    // disable UniProt headers until we know we have a hit
    for (std::string &output_path : mOutputPaths) {
        FS_dprint("DIAMOND file located at " + output_path + " being parsed");

        // reset uniprot info for each database
        is_uniprot = false;
        uniprot_attempts = 0;
        simSearchResults = {};
        ss.str("");
        ss.clear();

#ifdef TEST_SIM_PARSE_01
        test_ct = 0;
#endif

        // ensure file exists
        file_status = mpFileSystem->get_file_status(output_path);
        if (file_status & FileSystem::FILE_STATUS_EMPTY) {
            // Empty file, skip this file and let user know
            mpFileSystem->format_stat_stream(ss, "Compiled Similarity Search - DIAMOND - Best Overall");
            FS_dprint("WARNING: empty DIAMOND file found, skipping: " + output_path);
            ss << "DIAMOND file is empty. This will be skipped but pipeline will continue with other files";
            std::string out_msg = ss.str() + "\n";
            mpFileSystem->print_stats(out_msg);
            continue;   // WARNING CONTINUE if alignment file is empty

        } else if (file_status != 0) {
            throw ExceptionHandler("File not found or empty: " + output_path, ERR_ENTAP_RUN_SIM_SEARCH_FILTER);
        }

        // Begin using CSVReader lib to parse data
        io::CSVReader<DMND_COL_NUMBER, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(output_path);
        while (in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                           qstart, qend, sstart, send, evalue, bitscore, coverage,stitle)) {
            simSearchResults = {};

            // Get pointer to sequence in overall map
            QuerySequence *query = mpQueryData->get_sequence(qseqid);
            if (query == nullptr) {
                throw ExceptionHandler("Unable to find sequence in transcriptome: " + qseqid + " from file: " + output_path,
                                       ERR_ENTAP_RUN_SIM_SEARCH_FILTER);
            }

            // get species from database alignment title
            species = get_species(stitle);
            // get taxonomic information with species
            taxEntry = mpEntapDatabase->get_tax_entry(species);
            // get contaminant information
            contam_info = is_contaminant(taxEntry.lineage, mContaminateTaxons);

            // If this is a UniProt match and pull back info if so
            if (is_uniprot) {
                // Yes, UniProt match - Get uniprot info
                mpEntapDatabase->is_uniprot_entry(sseqid, simSearchResults.uniprot_info);
            } else {
                // No, not a UniProt database, check if it is
                if (uniprot_attempts <= UNIPROT_ATTEMPTS) {
                    // First UniProt match assumes the rest are UniProt as well in database
                    is_uniprot = mpEntapDatabase->is_uniprot_entry(sseqid, simSearchResults.uniprot_info);
                    if (!is_uniprot) {
                        uniprot_attempts++;
                    } else {
                        FS_dprint("Database file at " + output_path + "\nDetermined to be UniProt");
                        set_uniprot_headers();
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
            simSearchResults.is_informative = is_informative(stitle, mUninformativeTags);
            simSearchResults.is_informative ? simSearchResults.yes_no_inform = YES_FLAG :
                    simSearchResults.yes_no_inform  = NO_FLAG;

            query->add_alignment(mExecutionState, mSoftwareFlag,
                    simSearchResults, output_path, mInputLineage);

#ifdef TEST_SIM_PARSE_01
            if (++test_ct >= TEST_COUNT_MAX) {
                break;
            }
#endif
        } // END WHILE LOOP

        // Finished parsing and adding to alignment data, being to calc stats
        FS_dprint("File parsed, calculating statistics and writing output...");
        calculate_best_stats(false,output_path);
        FS_dprint("Success!");
    } // END FOR LOOP

    if (mParsedFile) {
        FS_dprint("Calculating overall Similarity Searching statistics...");
        calculate_best_stats(true);
        FS_dprint("Success!");
    } else {
        throw ExceptionHandler("No alignments found during Similarity Searching!",
                               ERR_ENTAP_RUN_SIM_SEARCH_FILTER);
    }
}

/**
 * ======================================================================
 * Function void ModDiamond::calculate_best_stats(bool is_final, std::string database_path)
 *
 * Description          - Calculates statistics from each database and a compiled
 *                        statistics
 *                      - Sends graphing information to manager to be graphed/printed
 *
 * Notes                - Empty database_path indicates overall database statistics
 *
 * @param is_final      - BOOL indicates overall statistics (TRUE) or individual database
 *                        statistics (FALSE)
 * @param database_path - Absolute path to database output to calculate statistics
 *                      - Empty path indicates "overall" statistics
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModDiamond::calculate_best_stats (bool is_final, std::string database_path) {

    GraphingManager::GraphingData graphingStruct;         // Graphing data
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
    uint64                      count_informative=0;    // Number of informative alignments
    uint64                      count_uninformative=0;  // Number of uninformative alignments
    uint64                      count_unselected=0;     // Number of unselected alignments (those that are not best hits)
    uint64                      count_TOTAL_alignments=0;
    uint32                      ct;
    fp64                        percent;
    fp64                        contam_percent;
    Compair<std::string>        contam_counter;
    Compair<std::string>        species_counter;
    Compair<std::string>        contam_species_counter;
    std::unordered_map<std::string, Compair<std::string>> frame_inform_counter;

    // Set up output directories for individual databases
    if (is_final) {
        // Overall results across databases
        base_path = mOverallResultsDir;
        database_shortname = "";
    } else {
        // Individual database results
        database_shortname = mPathToDatabase[database_path];
        base_path   = PATHS(mProcDir, database_shortname);
    }
    figure_base = PATHS(base_path, FIGURE_DIR);
    mpFileSystem->create_dir(base_path);
    mpFileSystem->create_dir(figure_base);

    // Open contam best hit tsv file and print headers
    std::string out_best_contams_filepath = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_HITS_CONTAM);
    mpQueryData->start_alignment_files(out_best_contams_filepath, mEntapHeaders, mGoLevels, mAlignmentFileTypes);

    // Open best hits files
    std::string out_best_hits_filepath = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_HITS);
    mpQueryData->start_alignment_files(out_best_hits_filepath, mEntapHeaders, mGoLevels, mAlignmentFileTypes);

    // Open best hits files with no contaminants
    std::string out_best_hits_no_contams = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_HITS_NO_CONTAM);
    mpQueryData->start_alignment_files(out_best_hits_no_contams, mEntapHeaders, mGoLevels, mAlignmentFileTypes);

    // Open unselected hits, so every hit that was not the best hit (tsv)
    std::string out_unselected_tsv  = PATHS(base_path, SIM_SEARCH_DATABASE_UNSELECTED);
    std::vector<FileSystem::ENT_FILE_TYPES> unselected_files = {FileSystem::ENT_FILE_DELIM_TSV};
    mpQueryData->start_alignment_files(out_unselected_tsv, mEntapHeaders, mGoLevels, unselected_files);

    // Open no hits file (fasta nucleotide)
    std::string out_no_hits_fa_nucl = PATHS(base_path, SIM_SEARCH_DATABASE_NO_HITS + FileSystem::EXT_FNN);
    std::ofstream file_no_hits_nucl(out_no_hits_fa_nucl, std::ios::out | std::ios::app);

    // Open no hits file (fasta protein)
    std::string out_no_hits_fa_prot  = PATHS(base_path, SIM_SEARCH_DATABASE_NO_HITS + FileSystem::EXT_FAA);
    std::ofstream file_no_hits_prot(out_no_hits_fa_prot, std::ios::out | std::ios::app);

    try {
        // Cycle through all sequences
        for (auto &pair : *mpQueryData->get_sequences_ptr()) {
            // Check if original sequences have hit a database
            if (!pair.second->hit_database(SIMILARITY_SEARCH, SIM_DIAMOND, database_path)) {
                // Did NOT hit a database during sim search
                // Do NOT log if it was never blasted
                if ((pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN) && mBlastp) ||
                    (!pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN) && !mBlastp)) {
                    // Protein/nucleotide did not hit database
                    count_no_hit++;
                    file_no_hits_nucl << pair.second->get_sequence_n() << std::endl;
                    file_no_hits_prot << pair.second->get_sequence_p() << std::endl;
                    // Graphing
                    frame = pair.second->getFrame();
                    if (!frame.empty()) {
                        if (frame_inform_counter.find(NO_HIT_FLAG) == frame_inform_counter.end()) {
                            frame_inform_counter.emplace(NO_HIT_FLAG, Compair<std::string>());
                        }
                        frame_inform_counter[NO_HIT_FLAG].add_value(frame);
                    }

                } else {
                    pair.second->set_blasted();
                }
            } else {
                // HIT a database during sim search

                QuerySequence::SimSearchResults *sim_search_data;
                SimSearchAlignment *best_hit;
                // Process unselected hits for non-final analysis and set best hit pointer
                if (is_final) {
                    best_hit =
                            pair.second->get_best_hit_alignment<SimSearchAlignment>(
                                    SIMILARITY_SEARCH, SIM_DIAMOND,"");
                    sim_search_data = best_hit->get_results();
                } else {
                    best_hit = pair.second->get_best_hit_alignment<SimSearchAlignment>(
                            SIMILARITY_SEARCH, SIM_DIAMOND,database_path);
                    QuerySequence::align_database_hits_t *alignment_data =
                            pair.second->get_database_hits(database_path,SIMILARITY_SEARCH, SIM_DIAMOND);
                    sim_search_data = best_hit->get_results();
                    for (auto &hit : *alignment_data) {
                        count_TOTAL_alignments++;
                        if (hit != best_hit) {  // If this hit is not the best hit
                            mpQueryData->add_alignment_data(out_unselected_tsv, pair.second, hit);
                            count_unselected++;
                        } else {
                            ;   // Do notthing
                        }
                    }
                }
                count_filtered++;   // increment best hit

                // Write to best hits files
                mpQueryData->add_alignment_data(out_best_hits_filepath, pair.second, best_hit);

                frame = pair.second->getFrame();     // Used for graphing
                species = sim_search_data->species;

                // Determine contaminant information and print to files
                if (sim_search_data->contaminant) {
                    // Species is considered a contaminant
                    count_contam++;
                    mpQueryData->add_alignment_data(out_best_contams_filepath, pair.second, best_hit);

                    contam = sim_search_data->contam_type;
                    contam_counter.add_value(contam);
                    contam_species_counter.add_value(species);
                } else {
                    // Species is NOT a contaminant, print to files
                    mpQueryData->add_alignment_data(out_best_hits_no_contams, pair.second, best_hit);
                }

                // Count species type
                species_counter.add_value(species);

                // Check if this is an informative alignment and respond accordingly
                if (sim_search_data->is_informative) {
                    count_informative++;
                    // Graphing
                    if (!frame.empty()) {
                        if (frame_inform_counter.find(INFORMATIVE_FLAG) == frame_inform_counter.end()) {
                            frame_inform_counter.emplace(INFORMATIVE_FLAG, Compair<std::string>());
                        }
                        frame_inform_counter[INFORMATIVE_FLAG].add_value(frame);
                    }

                } else {
                    count_uninformative++;
                    // Graphing
                    if (!frame.empty()) {
                        if (frame_inform_counter.find(UNINFORMATIVE_FLAG) == frame_inform_counter.end()) {
                            frame_inform_counter.emplace(UNINFORMATIVE_FLAG, Compair<std::string>());
                        }
                        frame_inform_counter[UNINFORMATIVE_FLAG].add_value(frame);
                    }
                }
            }
        }
    } catch (const std::exception &e){throw ExceptionHandler(e.what(), ERR_ENTAP_RUN_SIM_SEARCH_FILTER);}

    try {
        mpQueryData->end_alignment_files(out_best_contams_filepath);
        mpQueryData->end_alignment_files(out_best_hits_filepath);
        mpQueryData->end_alignment_files(out_best_hits_no_contams);
        mpQueryData->end_alignment_files(out_unselected_tsv);

        mpFileSystem->close_file(file_no_hits_nucl);
        mpFileSystem->close_file(file_no_hits_prot);
    } catch (const ExceptionHandler &e) {throw e;}

    // ------------ Calculate statistics and print to output ------------ //

    // ------------------- Setup graphing files ------------------------- //
    GraphingManager::GraphingData graph_contaminants_bar;
    GraphingManager::GraphingData graph_species_bar;
    graph_species_bar.x_axis_label   = "Species";
    graph_species_bar.y_axis_label   = "Count";
    graph_species_bar.text_file_path = PATHS(figure_base, GRAPH_SPECIES_BAR_TXT);
    graph_species_bar.fig_out_path   = PATHS(figure_base, GRAPH_SPECIES_BAR_PNG);
    graph_species_bar.graph_title    = database_shortname + GRAPH_SPECIES_TITLE;
    graph_species_bar.graph_type     = GraphingManager::ENT_GRAPH_BAR_HORIZONTAL;
    mpGraphingManager->initialize_graph_data(graph_species_bar);

    GraphingManager::GraphingData graph_frame_inform_stack;
    graph_frame_inform_stack.x_axis_label = "Category";
    graph_frame_inform_stack.y_axis_label = "Count";
    graph_frame_inform_stack.text_file_path = PATHS(figure_base, GRAPH_DATABASE_SUM_TXT);
    graph_frame_inform_stack.fig_out_path   = PATHS(figure_base, GRAPH_DATABASE_SUM_PNG);
    graph_frame_inform_stack.graph_title    = database_shortname + GRAPH_DATABASE_SUM_TITLE;
    graph_frame_inform_stack.graph_type     = GraphingManager::ENT_GRAPH_BAR_STACKED;
    mpGraphingManager->initialize_graph_data(graph_frame_inform_stack);


    // ------------------------------------------------------------------ //

    // Different headers if final analysis or database specific analysis
    if (is_final) {
        mpFileSystem->format_stat_stream(ss, "Compiled Similarity Search - DIAMOND - Best Overall");
    } else {
        mpFileSystem->format_stat_stream(ss, "Similarity Search - DIAMOND - " + database_shortname);
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
        mpFileSystem->print_stats(out_msg);
        return; // RETURN we do not have any alignments
    } else {
        mParsedFile = true;
    }

    // Sort counters
    contam_species_counter.sort(true);
    species_counter.sort(true);
    for (auto &pair : frame_inform_counter) {
        pair.second.sort(true);
    }

    ss <<
       "\n\tTotal unique transcripts with an alignment: " << count_filtered <<
       "\n\t\tReference transcriptome sequences with an alignment (FASTA):\n\t\t\t" << out_best_hits_filepath <<
       "\n\t\tSearch results (TSV):\n\t\t\t" << out_best_hits_filepath <<
       "\n\tTotal unique transcripts without an alignment: " << count_no_hit <<
       "\n\t\tReference transcriptome sequences without an alignment (FASTA):\n\t\t\t" << out_no_hits_fa_prot;
    // Have frame information
    if (frame_inform_counter.find(NO_HIT_FLAG) != frame_inform_counter.end()) {
        for (auto &pair : frame_inform_counter[NO_HIT_FLAG]._data) {
            ss << "\n\t\t" << pair.first << "(" << pair.second << ")";
            mpGraphingManager->add_datapoint(graph_frame_inform_stack.text_file_path, {pair.first, NO_HIT_FLAG,
                                                                                       std::to_string(pair.second)});
        }
    }
    ss <<
       "\n\tTotal unique informative alignments: " << count_informative;
    if (frame_inform_counter.find(INFORMATIVE_FLAG) != frame_inform_counter.end()) {
        for (auto &pair : frame_inform_counter[INFORMATIVE_FLAG]._data) {
            ss << "\n\t\t" << pair.first << "(" << pair.second << ")";
            mpGraphingManager->add_datapoint(graph_frame_inform_stack.text_file_path, {pair.first, INFORMATIVE_FLAG,
                                                                                       std::to_string(pair.second)});
        }
    }
    ss <<
       "\n\tTotal unique uninformative alignments: " << count_uninformative;
    if (frame_inform_counter.find(UNINFORMATIVE_FLAG) != frame_inform_counter.end()) {
        for (auto &pair : frame_inform_counter[UNINFORMATIVE_FLAG]._data) {
            ss << "\n\t\t" << pair.first << "(" << pair.second << ")";
            mpGraphingManager->add_datapoint(graph_frame_inform_stack.text_file_path, {pair.first, UNINFORMATIVE_FLAG,
                                                                                       std::to_string(pair.second)});
        }
    }

    // ********** Contaminant Calculations ************** //
    if (count_contam >= MIN_CONTAM_COUNT) {
        // Only show contaminant information if we have contaminants
        contam_percent = ((fp64) count_contam / count_filtered) * ENTAP_PERCENT;

        graph_contaminants_bar.x_axis_label   = "Contaminant Species";
        graph_contaminants_bar.y_axis_label   = "Count";
        graph_contaminants_bar.text_file_path = PATHS(figure_base, GRAPH_CONTAM_BAR_TXT);
        graph_contaminants_bar.fig_out_path   = PATHS(figure_base, GRAPH_CONTAM_BAR_PNG);
        graph_contaminants_bar.graph_title    = database_shortname + GRAPH_CONTAM_TITLE;
        graph_contaminants_bar.graph_type     = GraphingManager::ENT_GRAPH_BAR_HORIZONTAL;
        mpGraphingManager->initialize_graph_data(graph_contaminants_bar);

        ss <<
           "\n\tTotal unique contaminants: " << count_contam <<
           "(" << contam_percent << "%): " <<
           "\n\t\tTranscriptome reference sequences labeled as a contaminant (FASTA):\n\t\t\t"
           << out_best_contams_filepath <<
           "\n\t\tTranscriptome reference sequences labeled as a contaminant (TSV):\n\t\t\t" << out_best_contams_filepath;

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
            percent = ((fp64) pair.second / count_contam) * ENTAP_PERCENT;
            ss
                    << "\n\t\t\t" << ct << ")" << pair.first << ": "
                    << pair.second << "(" << percent << "%)";
            mpGraphingManager->add_datapoint(graph_contaminants_bar.text_file_path, {pair.first, std::to_string(pair.second)});
            ct++;
        }

        mpGraphingManager->graph_data(graph_contaminants_bar.text_file_path);
    }

    ss << "\n\tTop " << COUNT_TOP_SPECIES << " alignments by species:";
    ct = 1;
    for (auto &pair : species_counter._sorted) {
        if (ct > COUNT_TOP_SPECIES) break;
        percent = ((fp64) pair.second / count_filtered) * ENTAP_PERCENT;
        ss
                << "\n\t\t\t" << ct << ")" << pair.first << ": "
                << pair.second << "(" << percent << "%)";
        mpGraphingManager->add_datapoint(graph_species_bar.text_file_path, {pair.first, std::to_string(pair.second)});
        ct++;
    }
    std::string out_msg = ss.str() + "\n";
    mpFileSystem->print_stats(out_msg);
    // ------------------------------------------------------------------ //


    // -------------------------- Graphing Handle ----------------------- //
    mpGraphingManager->graph_data(graph_species_bar.text_file_path);
    mpGraphingManager->graph_data(graph_frame_inform_stack.text_file_path);
    // ------------------------------------------------------------------ //
}

void ModDiamond::set_uniprot_headers() {
    for (auto& header: UNIPROT_HEADERS) {
        mpQueryData->header_set(header, true);
    }
    mpQueryData->set_is_uniprot(true);
}

bool ModDiamond::set_version() {
    return false;
}
