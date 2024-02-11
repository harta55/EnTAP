/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2023, Alexander Hart, Dr. Jill Wegrzyn
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
#include "../ExceptionHandler.h"

std::vector<ENTAP_HEADERS> ModHorizontalGeneTransferDiamond::DEFAULT_HEADERS = {
        ENTAP_HEADER_HORIZONTALLY_TRANSFERRED_GENE
};

ModHorizontalGeneTransferDiamond::ModHorizontalGeneTransferDiamond(std::string &execution_stage_path, std::string &fasta_path,
                                                                   EntapDataPtrs &entap_data)
: AbstractHorizontalGeneTransfer(execution_stage_path, fasta_path, entap_data, "Horizontal Gene Transfer", DEFAULT_HEADERS){

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
            simSearchCmd.eval          = mEVal;
            simSearchCmd.tcoverage     = mTCoverage;
            simSearchCmd.qcoverage     = mQCoverage;
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

void ModHorizontalGeneTransferDiamond::parse() {

}

bool ModHorizontalGeneTransferDiamond::is_executable(std::string &exe) {
    return false;
}

bool ModHorizontalGeneTransferDiamond::set_version() {
    return false;
}

EntapModule::ModVerifyData ModHorizontalGeneTransferDiamond::verify_databases(ent_input_multi_str_t &databases, HGT_DATABASE_TYPES data_type) {
    EntapModule::ModVerifyData ret_data;
    std::string   database_name;        // Shortened name to be used for file naming
    std::string   out_path;             // Full output path for each database alignment
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

        // Check if file exists/can be read/empty
        file_status = mpFileSystem->get_file_status(out_path);
        if (file_status != 0) {
            FS_dprint("File for database " + database_name + " does not exist.\n" + out_path);
            // If we need to execute against ANY database
            ret_data.files_exist = false;
            // delete file just in case it is corrupt/empty
            mpFileSystem->delete_file(out_path);
        } else {
            // File found + is 'legit', can skip execution for it
            FS_dprint("File for database " + database_name + " exists, skipping...\n" + out_path);
            hgtDatabase.diamond_ran_success = true;
        }
        mHGTDatabases.push_back(hgtDatabase);
        ret_data.output_paths.push_back(out_path);   // Add paths to verify data (not currently used)
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
    return ret;}
