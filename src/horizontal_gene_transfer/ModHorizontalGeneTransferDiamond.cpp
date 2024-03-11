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

#include "ModHorizontalGeneTransferDiamond.h"
#include "csv.h"
#include "../QueryAlignment.h"

std::vector<ENTAP_HEADERS> ModHorizontalGeneTransferDiamond::DEFAULT_HEADERS = {
        ENTAP_HEADER_HORIZONTALLY_TRANSFERRED_GENE
};

ModHorizontalGeneTransferDiamond::ModHorizontalGeneTransferDiamond(std::string &execution_stage_path, std::string &fasta_path,
                                                                   EntapDataPtrs &entap_data)
: AbstractHorizontalGeneTransfer(execution_stage_path, fasta_path, entap_data, "HGT_DIAMOND", DEFAULT_HEADERS){

    mSoftwareFlag = HGT_DIAMOND;
    mExePath = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_DIAMOND_EXE);
}

EntapModule::ModVerifyData ModHorizontalGeneTransferDiamond::verify_files() {
    EntapModule::ModVerifyData donor_databases_verified = verify_databases(mDonorDatabasePaths, HGT_DATABASE_DONOR);
    EntapModule::ModVerifyData recipient_databases_verified = verify_databases(mRecipientDatabasePaths, HGT_DATABASE_RECIPIENT);
    EntapModule::ModVerifyData return_databases_verified;

    return_databases_verified.files_exist = donor_databases_verified.files_exist && recipient_databases_verified.files_exist;

    return return_databases_verified;
}

/***
 * Horizontal Gene Transfer
 * This is a process of determining if a gene has been horizontally transferred.
 *  There are several stages to be performed:
 *
 *  1. Run Similarity Search with our transcriptome against all of our Donor and Recipient Databases
 *  2. Leverage alignment information to determine which genes are 'HGT candidates'. These 'may' be HGT genes, but
 *      we'll need to use more information to be sure. This is where the GFF file comes into play
 *  3. Using GFF input from user, determine upstream and downstream genes to our 'HGT candidates'
 *  4. If upstream/downstream genes match the taxonomy of our Candidate then this is a true HGT gene
 */
void ModHorizontalGeneTransferDiamond::execute() {
    uint16 file_status = 0;     // File statuses (empty, cannot be read, etc...)
    AbstractSimilaritySearch::SimSearchCmd simSearchCmd;  // DIAMOND commands

    // 1. Similarity Search against Donor and Recipient databases
    FS_dprint("HGT Running Similarity Search against Donor and Recipient Databases");
    for (auto &database : mHGTDatabases) {
        file_status = mpFileSystem->get_file_status(database.diamond_output);
        if (file_status != 0) {
            // If file does not exist or cannot be read, execute diamond
            FS_dprint("File not found, executing against database at: " + database.database_path);

            simSearchCmd = {};
            simSearchCmd.database_path = database.database_path;
            simSearchCmd.output_path   = database.diamond_output;
            simSearchCmd.std_out_path  = database.diamond_output + FileSystem::EXT_STD;
            simSearchCmd.threads       = (uint16)mThreads;
            simSearchCmd.query_path    = mInputTranscriptome;
            simSearchCmd.eval          = DMND_E_VAL;
            simSearchCmd.tcoverage     = DMND_TARGET_COVERAGE;
            simSearchCmd.qcoverage     = DMND_QUERY_COVERAGE;
            simSearchCmd.exe_path      = mExePath;
            simSearchCmd.blastp        = mBlastp;

            try {
                run_diamond(&simSearchCmd);
                database.diamond_ran_success = true;
            } catch (const ExceptionHandler &e ){
                throw e;
            }
            FS_dprint("Success! Results written to: " + database.diamond_output);
        }
    }
    FS_dprint("HGT run against donor and recipients complete, determining HGT candidates...");
    // 2.
}

// TODO restrucutre this and ModDIAMOND so that code can be reused better
void ModHorizontalGeneTransferDiamond::parse() {
    uint16              file_status=0;                  // File statuses (defined FileSystem.h)
    std::string         species;                        // Species from subject sequence
    QuerySequence::HorizontalGeneTransferResults horizontalGeneTransferResults;   // Compiled similarity search results
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

    FS_dprint("Beginning to filter individual DIAMOND files...");
    TC_print(TC_PRINT_COUT, "Parsing DIAMOND Horizontal Gene Transfer...");

    for (HGTDatabase &hgtDatabase : mHGTDatabases) {
        FS_dprint("DIAMOND file located at " + hgtDatabase.diamond_output + " being parsed");

        horizontalGeneTransferResults = {};
        ss.str("");
        ss.clear();

        // ensure file exists
        file_status = mpFileSystem->get_file_status(hgtDatabase.diamond_output);
        if (file_status != 0) {
            // Empty file, skip this file and let user know
            mpFileSystem->format_stat_stream(ss, "Compiled Horizontal Gene Transfer - DIAMOND - Best Overall");
            FS_dprint("WARNING: empty DIAMOND file found, skipping: " + hgtDatabase.diamond_output);
            ss << "DIAMOND file is empty or cannot be read. This will be skipped but pipeline will continue with other files";
            std::string out_msg = ss.str() + "\n";
            mpFileSystem->print_stats(out_msg);
            continue;   // WARNING CONTINUE if alignment file is empty
        }

        // Begin using CSVReader lib to parse data
        io::CSVReader<DMND_COL_NUMBER, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(hgtDatabase.diamond_output);
        while (in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                           qstart, qend, sstart, send, evalue, bitscore, coverage,stitle)) {
            horizontalGeneTransferResults = {};

            // Get pointer to sequence in overall map
            QuerySequence *query = mpQueryData->get_sequence(qseqid);
            if (query == nullptr) {
                throw ExceptionHandler("Unable to find sequence in transcriptome (" + qseqid + ") from file: " + hgtDatabase.diamond_output,
                                       ERR_ENTAP_HGT_PARSE);
            }

            // get species from database alignment title
            species = AbstractSimilaritySearch::get_species(stitle);
            // get taxonomic information with species
            taxEntry = mpEntapDatabase->get_tax_entry(species);

            // Compile sim search data
            horizontalGeneTransferResults.database_path = hgtDatabase.database_path;
            horizontalGeneTransferResults.qseqid = qseqid;
            horizontalGeneTransferResults.sseqid = sseqid;
            horizontalGeneTransferResults.pident = pident;
            horizontalGeneTransferResults.length = length;
            horizontalGeneTransferResults.mismatch = mismatch;
            horizontalGeneTransferResults.gapopen = gapopen;
            horizontalGeneTransferResults.qstart = qstart;
            horizontalGeneTransferResults.qend = qend;
            horizontalGeneTransferResults.sstart = sstart;
            horizontalGeneTransferResults.send = send;
            horizontalGeneTransferResults.stitle = stitle;
            horizontalGeneTransferResults.bit_score = bitscore;
            horizontalGeneTransferResults.lineage = taxEntry.lineage;
            horizontalGeneTransferResults.species = species;
            horizontalGeneTransferResults.e_val_raw = evalue;
            horizontalGeneTransferResults.e_val = float_to_sci(evalue,2);
            horizontalGeneTransferResults.coverage_raw = coverage;
            horizontalGeneTransferResults.coverage = float_to_string(coverage);

            query->add_alignment(mExecutionState, mSoftwareFlag,
                                 horizontalGeneTransferResults, hgtDatabase.database_path);

        } // END WHILE LOOP

        // Finished parsing DIAMOND output file and adding to alignment data, begin to calc stats
        FS_dprint("File parsed, calculating statistics and writing output...");
        calculate_best_stats(hgtDatabase);
        FS_dprint("Success!");
    } // END FOR LOOP
    calculate_hgt_candidates(mHGTDatabases);

    TC_print(TC_PRINT_COUT, "Success");
}

bool ModHorizontalGeneTransferDiamond::is_executable(std::string &exe) {
    TerminalData terminalData;  // Terminal data structure

    terminalData.command = exe + " --version";
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;

    return TC_execute_cmd(terminalData) == 0;
}

bool ModHorizontalGeneTransferDiamond::set_version() {
    return false;
}

EntapModule::ModVerifyData ModHorizontalGeneTransferDiamond::verify_databases(ent_input_multi_str_t &databases, HGT_DATABASE_TYPES data_type) {
    EntapModule::ModVerifyData ret_data;
    uint16 file_status = 0;             // File statuses (empty, doesn't exist, etc...)

    ret_data.files_exist = true;

    for (std::string &data_path : databases) {
        FS_dprint("Verifying previous execution of database: " + data_path + "...");

        HGTDatabase hgtDatabase;

        hgtDatabase.database_shortname = get_database_shortname(data_path);
        hgtDatabase.database_path = data_path;
        hgtDatabase.diamond_output = get_database_output_path(data_path);
        hgtDatabase.diamond_ran_success = false;
        hgtDatabase.database_type = data_type;

        if (data_type == HGT_DATABASE_RECIPIENT) mRecipientDatabaseCt++;
        if (data_type == HGT_DATABASE_DONOR) mDonorDatabaseCt++;

        // Check if file exists/can be read/empty
        file_status = mpFileSystem->get_file_status(hgtDatabase.diamond_output);
        if (file_status != 0) {
            FS_dprint("File for database " + hgtDatabase.database_shortname + " does not exist at.\n" + hgtDatabase.diamond_output);
            // If we need to execute against ANY database
            ret_data.files_exist = false;
            // delete file just in case it is corrupt/empty
            mpFileSystem->delete_file(hgtDatabase.diamond_output);
        } else {
            // File found + is 'legit', can skip execution for it
            FS_dprint("File for database " + hgtDatabase.database_shortname + " exists, skipping...\n" + hgtDatabase.diamond_output);
            hgtDatabase.diamond_ran_success = true;
        }
        mHGTDatabases.push_back(hgtDatabase);
        ret_data.output_paths.push_back(hgtDatabase.diamond_output);   // Add paths to verify data (not currently used)
    }
    FS_dprint("Success! Verified files for DIAMOND, continuing...");
    return ret_data;
}

bool ModHorizontalGeneTransferDiamond::run_diamond(AbstractSimilaritySearch::SimSearchCmd *cmd) {
    TerminalData    terminalData;   // Terminal data
    command_map_t   tc_commands;    // Terminal command map
    int32           err_code;       // Error codee from terminal execution
    bool            ret = true;     // Return value, if execution has succeeded
    std::string     temp_exe;       // DIAMOND needs blastp/x directly after DIAMOND exe

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

void ModHorizontalGeneTransferDiamond::calculate_best_stats (ModHorizontalGeneTransferDiamond::HGTDatabase &database) {

    GraphingManager::GraphingData graphingStruct;         // Graphing data
    std::string                 species;
    std::string                 frame;
    std::string                 base_path;
    std::string                 temp_file_path;
    std::string                 contam;
    std::stringstream           ss;
    uint64                      count_no_hit=0;
    uint64                      count_filtered=0;
    uint64                      count_unselected=0;     // Number of unselected alignments (those that are not best hits)
    uint64                      count_TOTAL_alignments=0;
    uint32                      ct;
    fp64                        percent;
    Compair<std::string>        species_counter;
    std::unordered_map<std::string, Compair<std::string>> frame_inform_counter;

    base_path   = PATHS(mProcDir, database.database_shortname);
    mpFileSystem->create_dir(base_path);

    // Open best hits files
    std::string out_best_hits_filepath = PATHS(base_path, DIAMOND_PREFIX + HGT_DIAMOND_DATABASE_BEST_HITS);
    mpQueryData->start_alignment_files(out_best_hits_filepath, mEntapHeaders, mGoLevels, mAlignmentFileTypes);

    // Open unselected hits, so every hit that was not the best hit (tsv)
    std::string out_unselected_tsv  = PATHS(base_path, DIAMOND_PREFIX + HGT_DIAMOND_DATABASE_UNSELECTED);
    std::vector<FileSystem::ENT_FILE_TYPES> unselected_files = {FileSystem::ENT_FILE_DELIM_TSV};
    mpQueryData->start_alignment_files(out_unselected_tsv, mEntapHeaders, mGoLevels, unselected_files);

    // No hits files - only print what we have data for
    std::vector<FileSystem::ENT_FILE_TYPES> no_hits_files;
    std::string out_no_hits = PATHS(base_path, HGT_DIAMOND_DATABASE_NO_HITS);
    if (mpQueryData->DATA_FLAG_GET(QueryData::IS_PROTEIN)) {
        no_hits_files.push_back(FileSystem::ENT_FILE_FASTA_FAA);
    }
    if (mpQueryData->DATA_FLAG_GET(QueryData::IS_NUCLEOTIDE)) {
        no_hits_files.push_back(FileSystem::ENT_FILE_FASTA_FNN);
    }
    if (!no_hits_files.empty()) {
        mpQueryData->start_alignment_files(out_no_hits, mEntapHeaders, mGoLevels, no_hits_files);
    }

    try {
        // Cycle through all sequences
        for (auto &pair : *mpQueryData->get_sequences_ptr()) {
            // Check if original sequences have hit a database
            if (!pair.second->hit_database(HORIZONTAL_GENE_TRANSFER, HGT_DIAMOND, database.database_path)) {
                // Did NOT hit a database during sim search
                // Do NOT log if it was never blasted
                if ((pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN) && mBlastp) ||
                    (!pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN) && !mBlastp)) {
                    // Protein/nucleotide did not hit database
                    count_no_hit++;
                    mpQueryData->add_alignment_data(out_no_hits, pair.second, nullptr); // Function checks whether files been initializd

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
                // HIT a database during HGT DIAMOND search

                QuerySequence::HorizontalGeneTransferResults *hgt_data;
                HorizontalGeneTransferDmndAlignment *best_hit;
                // Process unselected hits for non-final analysis and set best hit pointer

                best_hit = pair.second->get_best_hit_alignment<HorizontalGeneTransferDmndAlignment>(
                        HORIZONTAL_GENE_TRANSFER, HGT_DIAMOND,database.database_path);
                QuerySequence::align_database_hits_t *alignment_data =
                        pair.second->get_database_hits(database.database_path,HORIZONTAL_GENE_TRANSFER, HGT_DIAMOND);
                hgt_data = best_hit->get_results();
                for (auto &hit : *alignment_data) {
                    count_TOTAL_alignments++;
                    if (hit != best_hit) {  // If this hit is not the best hit
                        mpQueryData->add_alignment_data(out_unselected_tsv, pair.second, hit);
                        count_unselected++;
                    } else {
                        ;   // Do notthing
                    }
                }

                count_filtered++;   // increment best hit

                // Write to best hits files
                mpQueryData->add_alignment_data(out_best_hits_filepath, pair.second, best_hit);

                frame = pair.second->getFrame();     // Used for graphing
                species = hgt_data->species;

                // Count species type
                species_counter.add_value(species);
            }
        }
    } catch (const std::exception &e){throw ExceptionHandler(e.what(), ERR_ENTAP_HGT_PARSE);}

    try {
        mpQueryData->end_alignment_files(out_best_hits_filepath);
        mpQueryData->end_alignment_files(out_unselected_tsv);
        mpQueryData->end_alignment_files(out_no_hits);
    } catch (const ExceptionHandler &e) {throw e;}

    // ------------ Calculate statistics and print to output ------------ //
    // ------------------------------------------------------------------ //

    // Different headers if final analysis or database specific analysis

    mpFileSystem->format_stat_stream(ss, "Horizontal Gene Transfer - DIAMOND - " + database.database_shortname);
    ss <<
       "Search results:\n"                    << database.database_path <<
       "\n\tTotal alignments: "               << count_TOTAL_alignments   <<
       "\n\tTotal unselected results: "       << count_unselected      <<
       "\n\t\tWritten to: "                   << out_unselected_tsv;

    // If no total or file alignments for this database, return and warn user
    if (count_TOTAL_alignments == 0 || count_filtered == 0) {
        ss << "WARNING: No alignments for this database";
        std::string out_msg = ss.str() + "\n";
        mpFileSystem->print_stats(out_msg);
        return; // RETURN we do not have any alignments
    }

    // Sort counters
    species_counter.sort(true);
    for (auto &pair : frame_inform_counter) {
        pair.second.sort(true);
    }

    ss <<
       "\n\tTotal unique transcripts with an alignment: " << count_filtered <<
       "\n\t\tReference transcriptome sequences with an alignment (FASTA):\n\t\t\t" << out_best_hits_filepath <<
       "\n\t\tSearch results (TSV):\n\t\t\t" << out_best_hits_filepath <<
       "\n\tTotal unique transcripts without an alignment: " << count_no_hit <<
       "\n\t\tReference transcriptome sequences without an alignment (FASTA):\n\t\t\t" << out_no_hits;
    // Have frame information
    if (frame_inform_counter.find(NO_HIT_FLAG) != frame_inform_counter.end()) {
        for (auto &pair : frame_inform_counter[NO_HIT_FLAG]._data) {
            ss << "\n\t\t" << pair.first << "(" << pair.second << ")";
        }
    }

    ss << "\n\tTop " << COUNT_TOP_SPECIES << " alignments by species:";
    ct = 1;
    for (auto &pair : species_counter._sorted) {
        if (ct > COUNT_TOP_SPECIES) break;
        percent = ((fp64) pair.second / count_filtered) * ENTAP_PERCENT;
        ss
                << "\n\t\t\t" << ct << ")" << pair.first << ": "
                << pair.second << "(" << percent << "%)";
        ct++;
    }
    std::string out_msg = ss.str() + "\n";
    mpFileSystem->print_stats(out_msg);
}

void ModHorizontalGeneTransferDiamond::calculate_hgt_candidates(std::vector<HGTDatabase> hgt_databases) {
    FS_dprint("Calculating HGT Candidates...");

    uint16 query_donor_ct=0;
    uint16 query_recipient_ct=0;
    const QuerySequence *upstream_sequence;
    const QuerySequence *downstream_sequence;

    /*
     * 1) Determine Horizontal Gene Transfer Candidates
     *
     *  For this stage, we cycle through the entire transcriptome to determine
     *      the HGT candidates based on the Donor and Recipient alignments
     *
     *  HGT Candidate Determined By:
     *      1. One or more hits against the DONOR database
     *      2. No hits against RECIPIENT databases
     *
     * */
    // Loop through entire transcriptome, probably slow TODO speed up HGT parsing
    for (auto &pair : *mpQueryData->get_sequences_ptr()) {

        // Check if sequence hit at least one HGT database
        if (pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_HGT_BLASTED)) {
            query_donor_ct=0;
            query_recipient_ct=0;
            for (auto &database : hgt_databases) {
                // Check if this query hit against the database then keep track of count
                if (pair.second->hit_database(HORIZONTAL_GENE_TRANSFER, HGT_DIAMOND,database.database_path)) {
                    if (database.database_type == HGT_DATABASE_DONOR) {
                        query_donor_ct++;
                    } else {
                        query_recipient_ct++;
                    }
                }
            }
            pair.second->setMDonorDatabaseHitCt(query_donor_ct);
            pair.second->setMRecipientDatabaseHitCt(query_recipient_ct);

            // Determine if horizontal gene transfer candidate
            // If 1 or more hits against DONOR databases BUT not all
            if ((query_donor_ct >= HGT_DONOR_DATABASE_MIN) && (query_donor_ct < mDonorDatabaseCt)) {
                // AND we have no recipient hits
                if (query_recipient_ct == 0) {
                    pair.second->QUERY_FLAG_CHANGE(QuerySequence::QUERY_HGT_CANDIDATE, true);
                    FS_dprint("HGT Candidate Found! " + pair.first);
                }
            }
        }
    }

    /*
     * 2) Refine HGT Candidates based on neighboring genes
     *
     *  For this stage, we cycle through our HGT Candidate and analyze relationships between
     *      the neighbors (using GFF file information)
     *
     *  Refinement determined by:
     *      1. If either neighbor of an HGT Candidate is an HGT Candidate, remove them all from consideration
     *      2. Neighbors should not have any DONOR database hits, only RECIPIENT
     * */
    for (auto &pair : *mpQueryData->get_sequences_ptr()) {

        // We only want to analyze those sequences that are an HGT Candidate
        if (pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_HGT_CANDIDATE)) {

            // 2.1) Confirm we have existing neighbors, if not, NOT HGT
            upstream_sequence = pair.second->getMpUpstreamSequence();
            downstream_sequence = pair.second->getMpDownstreamSequence();
            // Make sure we have valid upstream and downstream genes first
            if ((upstream_sequence == nullptr) || (downstream_sequence == nullptr)) {
                FS_dprint("WARNING unable to verify gene (" + pair.first + ") since upstream/downstream gene missing");
                continue;   // !!!! WARNING CONTINUE, gene neighbors missing
            }

            // 2.2) Check if neighbors are HGT Candidates. If so, NOT HGT
            if ((upstream_sequence->QUERY_FLAG_GET(QuerySequence::QUERY_HGT_CANDIDATE)) ||
                (downstream_sequence->QUERY_FLAG_GET(QuerySequence::QUERY_HGT_CANDIDATE))) {

                pair.second->QUERY_FLAG_CHANGE(QuerySequence::QUERY_HGT_CONFIRMED, false);
                FS_dprint("WARNING skipped gene (" + pair.first + ") due to neighbor HGT Candidates");
                continue; // !!!! WARNING CONTINUE, neighbor HGT Candidates
            }

            // 2.3) Check if neighbors have any donor hits, if so, NOT HGT
            //  neighbors should ONLY have recipient hits
            if ((upstream_sequence->getMDonorDatabaseHitCt() >= HGT_DONOR_DATABASE_NEIGHBOR_MAX) ||
                (downstream_sequence->getMDonorDatabaseHitCt() >= HGT_DONOR_DATABASE_NEIGHBOR_MAX)) {

                pair.second->QUERY_FLAG_CHANGE(QuerySequence::QUERY_HGT_CONFIRMED, false);
                FS_dprint("WARNING skipped gene (" + pair.first + ") due to neighbor HGT Candidates having donor hits");
                continue; // !!!! WARNING CONTINUE, neighbor HGT Candidates hit donor database
            }

            // If we got past every other step, this is HGT Confirmed!!
            pair.second->QUERY_FLAG_CHANGE(QuerySequence::QUERY_HGT_CONFIRMED, true);
            FS_dprint("HGT CANDIDATE CONFIRMED!!!: " + pair.first);
        }
    }

    // Loop through only our HGT Candidates for further analyze them

    FS_dprint("HGT Candidate calculation complete!");
}
