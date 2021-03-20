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
 * MERCHANTABILIT * but WITHOUT ANY WARRANTY; without even the implied warranty of
Y or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/


//*********************** Includes *****************************
#include "../common.h"
#include "ModGeneMarkST.h"
#include "../ExceptionHandler.h"
#include "../FileSystem.h"
#include "../TerminalCommands.h"
#include "../GraphingManager.h"
//**************************************************************

std::vector<ENTAP_HEADERS> ModGeneMarkST::DEFAULT_HEADERS = {
        ENTAP_HEADER_FRAME
};


/**
 * ======================================================================
 * Function std::pair<bool, std::string> ModGeneMarkST::verify_files()
 *
 * Description           - Checks whether GeneMarkS-T has already been ran
 *                         with the same input
 *
 * Notes                 - None
 *
 *
 * @return               - Bool   - Yes if previous files have been found
 *                       - String - Path to frame selected transcriptome
 *
 * =====================================================================
 */
EntapModule::ModVerifyData ModGeneMarkST::verify_files() {
    ModVerifyData modVerifyData;    // Data to return from verification process

    modVerifyData.files_exist = false;

    FS_dprint("Beginning to verify GeneMarkS-T module files...");

    // If GeneMarkS-T has already been ran
    if (mpFileSystem->file_exists(mFinalFaaPath) && mpFileSystem->file_exists(mFinalLstPath)) {
        // Yes, files have been found and GeneMarkS-T has already been ran
        FS_dprint("Files found at: " + mFinalFaaPath + "\nand: " + mFinalLstPath +
                "\ncontinuing EnTAP with these files and skipping Frame Selection");
        modVerifyData.files_exist = true;
    } else {
        // No, files not found. GeneMarkS-T needs to be ran
        FS_dprint("File not found at " + mFinalFaaPath + "\nor " + mFinalLstPath +
                  " so continuing with Frame Selection");
    }
    modVerifyData.output_paths = vect_str_t{mFinalFaaPath};
    return modVerifyData;
}


/**
 * ======================================================================
 * Function std::string ModGeneMarkST::execute()
 *
 * Description          - Will run GeneMarkS-T through PStreams terminal
 *                        commands
 *                      - 4/5 output files are expected (3 are required to continue to parsing)
 *
 * Notes                - Fatal Error if an issue occurred during GeneMarkS-T execution
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModGeneMarkST::execute() {
    // Outfiles: file/path.faa, file/path.fnn
    // WARNING GeneMarkS-T assumes working directory as output right now
    std::string     temp_lst_file;      // Absolute path to output lst file from GeneMark (cwd)
    std::string     temp_faa_file;      // Absolute path to output FAA file from GeneMark (cwd)
    std::string     temp_fnn_file;      // Absolute path to output FNN file from GeneMark (cwd)
    std::string     temp_hmm_file;      // Absolute path to output HMM file from GeneMark (cwd)
    std::string     temp_gms_log_file;  // Absolute path to output log file from GeneMark (cwd)
    std::string     genemark_cmd;       // Terminal comand for GeneMarkS-T
    std::string     genemark_std_out;   // Standard output from GeneMark terminal run
    std::string     line;               // Temp line string
    int32           err_code;           // Terminal command error code
    TerminalData    terminalData;       // Terminal command structure (Defined in TerminalCommands.h)

    // Set temporary outputs (these will be moved to FINAL outpaths
    temp_faa_file     = PATHS(mpFileSystem->get_cur_dir(), mTranscriptomeFilename + FileSystem::EXT_FAA);
    temp_fnn_file     = PATHS(mpFileSystem->get_cur_dir(), mTranscriptomeFilename + FileSystem::EXT_FNN);
    temp_lst_file     = PATHS(mpFileSystem->get_cur_dir(), mTranscriptomeFilename + FileSystem::EXT_LST);
    temp_hmm_file     = PATHS(mpFileSystem->get_cur_dir(), GENEMARK_HMM_FILE);
    temp_gms_log_file = PATHS(mpFileSystem->get_cur_dir(), GENEMARK_LOG_FILE);

    genemark_cmd     = mExePath + " -faa -fnn " + mInputTranscriptome;
    genemark_std_out = PATHS(mModOutDir, GENEMARK_STD_OUT);

    terminalData.command        = genemark_cmd;
    terminalData.print_files    = true;
    terminalData.suppress_std_err = false;
    terminalData.base_std_path  = genemark_std_out;

    // !!!WARNING!!! GeneMarkS-T output is always in the CWD
    err_code = TC_execute_cmd(terminalData);
    if (err_code != 0 ) {
        throw ExceptionHandler("Error in running GeneMarkST at file located at: " +
                               mInputTranscriptome + "\nGeneMarkST Error:\n" + terminalData.err_stream,
                               ERR_ENTAP_RUN_GENEMARK);
    }
    FS_dprint("Success!");

    // Ensure files successfully printed
    if (!mpFileSystem->file_exists(temp_faa_file)) {
        throw ExceptionHandler("Error unable to find GeneMarkS-T file located at: " + temp_faa_file,
                               ERR_ENTAP_RUN_GENEMARK);
    } else if (!mpFileSystem->file_exists(temp_fnn_file)) {
        throw ExceptionHandler("Error unable to find GeneMarkS-T file located at: " + temp_fnn_file,
                               ERR_ENTAP_RUN_GENEMARK);
    } else if (!mpFileSystem->file_exists(temp_lst_file)) {
        throw ExceptionHandler("Error unable to find GeneMarkS-T file located at: " + temp_lst_file,
                               ERR_ENTAP_RUN_GENEMARK);
    }

    FS_dprint("GeneMarkS-T files printed to:\nFAA File: " + temp_faa_file +
                                            "\nFNN File: " + temp_fnn_file+
                                            "\nLST File: " + temp_lst_file);

    // No longer doing this, just changing CWD
#if 0
    // Format genemarks-t output (remove blank lines for user)
    // Format FNN file
    FS_dprint("Formatting GeneMarkST files...");
    try {
        std::ifstream in_file_fnn(temp_fnn_file);
        std::ofstream out_file_fnn(mFinalFnnPath);
        while(getline(in_file_fnn, line)) {
            if (!line.empty()) {
                out_file_fnn << line << '\n';
            }
        }
        in_file_fnn.close();
        out_file_fnn.close();

        // Format FAA file
        std::ifstream in_file_faa(temp_faa_file);
        std::ofstream out_file_faa(mFinalFaaPath);
        while(getline(in_file_faa, line)) {
            if (!line.empty()) {
                out_file_faa << line << '\n';
            }
        }
        in_file_faa.close();
        out_file_faa.close();
    } catch (std::exception &e) {
        throw ExceptionHandler("Error formatting GeneMarkS-T results\n" + std::string(e.what()),
                               ERR_ENTAP_RUN_GENEMARK_MOVE);
    }

    // Delete temporary files (ignore errors)
    mpFileSystem->delete_file(temp_fnn_file);
    mpFileSystem->delete_file(temp_faa_file);
#endif

    // Move other output files to module directory
    if (!mpFileSystem->rename_file(temp_lst_file, mFinalLstPath)) {
        throw ExceptionHandler("Error moving GeneMarkS-T results", ERR_ENTAP_RUN_GENEMARK_MOVE);
    }

    // move log file, ignore errors not needed for execution
    mpFileSystem->rename_file(temp_gms_log_file, mFinalGmstLogPath);

    // Does not always exist, not needed for final calculation so ignore errors
    mpFileSystem->rename_file(temp_hmm_file,mFinalHmmPath);
    FS_dprint("Success!");
}


/**
 * ======================================================================
 * Function void ModGeneMarkST::parse()
 *
 * Description          - Parses/calculates statistics of GeneMark output
 *                        and prints for user
 *                      - Calls Graphing Manager to graph results
 *                      - Prints genes (partial, complete, internal...) to different files
 *                      - Updates QueryData structure with new frame information (sequence, frame)
 *
 * Notes                - Fatal Error if any issue occurs
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModGeneMarkST::parse() {
    // generate maps, query->sequence
    FS_dprint("Beginning to calculate Genemark statistics...");

    // Ensure paths we need exist
    if (!mpFileSystem->file_exists(mFinalFaaPath)) {
        throw ExceptionHandler("Final GeneMarkST output not found at: " + mFinalFaaPath,
            ERR_ENTAP_RUN_GENEMARK_PARSE);
    } else if (!mpFileSystem->file_exists(mFinalLstPath)) {
        throw ExceptionHandler("Final GeneMarkST lst output not found at: " + mFinalLstPath,
                               ERR_ENTAP_RUN_GENEMARK_PARSE);
    } else if (!mpFileSystem->file_exists(mFinalFnnPath)) {
        throw ExceptionHandler("Final GeneMarkST lst output not found at: " + mFinalFnnPath,
                               ERR_ENTAP_RUN_GENEMARK_PARSE);
    }

    try {
        // ------------------- Parse GeneMarkS-T Files ---------------------- //
        FS_dprint("Beginning to parse GeneMarkS-T files...");
        // Parse protein file (.faa)
        genemark_parse_fasta(mFinalFaaPath, FileSystem::ENT_FILE_FASTA_FAA);
        // Parse nucleotide file (.fnn) TODO cleanup/ move to query data just throwing in for now
        genemark_parse_fasta(mFinalFnnPath, FileSystem::ENT_FILE_FASTA_FNN);
        // Parse lst file to get info for each sequence (partial, internal...)
        genemark_parse_lst(mFinalLstPath);
        FS_dprint("Success! GeneMarkS-T files parsed");
        // ------------------------------------------------------------------ //

        frame_calculate_statistics();

    } catch (const ExceptionHandler &e) {throw e;}
    FS_dprint("Success! Parsing complete");
}


/**
 * ======================================================================
 * Function void ModGeneMarkST::genemark_parse_protein
 *                                                  (std::string &protein)
 *
 * Description          - Analyzes fasta file to pull protein/nucleotide sequences that
 *                        will be added to QueryData
 *
 * Notes                - Updates QueryData structure with nucleotide or protien
 *
 * @param protein       - Absolute path to fasta (.faa or .fnn) file
 * @param file_type     - Type of file to parse, either faa or fnn
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModGeneMarkST::genemark_parse_fasta(std::string &fasta, FileSystem::ENT_FILE_TYPES file_type) {
    FS_dprint("Parsing GeneMarkS-T FASTA file at: " + fasta);

    std::string     line;               // Temp string
    std::string     sequence;           // FAA or FNN of translated sequence
    std::string     seq_id;             // ID of translated sequence
    uint16          first;
    uint16          second;
    QuerySequence  *pQuery_sequence;

    // Filepath checked before
    std::ifstream in_file(fasta);
    while(true) {
        getline(in_file,line);
        if (line.empty() && !in_file.eof()) continue;
        if (line.find(FileSystem::FASTA_FLAG)==0 || in_file.eof()) {
            if (!seq_id.empty()) {
                if (in_file.eof()) {
                    sequence += line + "\n";
                }
                pQuery_sequence = mpQueryData->get_sequence(seq_id);
                if (pQuery_sequence == nullptr) {
                    throw ExceptionHandler("Unable to find correct sequence ID in GeneMarkS-T run", ERR_ENTAP_RUN_GENEMARK_PARSE);
                } else {
                    if (file_type == FileSystem::ENT_FILE_FASTA_FAA) {
                        pQuery_sequence->set_sequence_p(sequence);
                    } else if (file_type == FileSystem::ENT_FILE_FASTA_FNN) {
                        pQuery_sequence->set_sequence_n(sequence);
                    }
                }
            }
            if (in_file.eof()) break;
            // GeneMarkS-T adds a tab character within the sequence headers
            // Remove this to find the 'actual' sequence header we can reference against
            // our original transcriptome
            first    = (uint16) (line.find(FileSystem::FASTA_FLAG)+1);
            second   = (uint16) line.find('\t');
            if (first != std::string::npos || second != std::string::npos) {
                seq_id   = line.substr(first,line.find('\t')-first);
                sequence = FileSystem::FASTA_FLAG + seq_id + "\n";
            }
        } else {
            sequence += line + "\n";
        }
    }

    FS_dprint("Success! File parsed");
    in_file.close();
}


/**
 * ======================================================================
 * Function void ModGeneMarkST::genemark_parse_lst(std::string &lst_path)
 *
 * Description          - Parses .lst file produced from GeneMark run to get Frame types
 *                        (partial, complete, etc...)
 *                      - Will update QueryData with frame information
 *
 * Notes                - Updates mpQueryData
 *
 * @param lst_path     - Path to .lst file produced from GeneMark
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModGeneMarkST::genemark_parse_lst(std::string &lst_path) {
    FS_dprint("Parsing GeneMarkS-T LST file at: " + lst_path);

    std::string     frame;
    std::string     line;
    std::string     seq_id;
    bool            prime_5;
    bool            prime_3;
    uint16          first;
    QuerySequence   *pQuery_sequence;

    // Filepath checked before
    std::ifstream in_file(lst_path);
    while (getline(in_file,line)) {
        if (line.empty()) continue;
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        if (line.find("FASTA") == 0) {
            first = (uint16) (line.find(':') + 1);
            seq_id = line.substr(first);
        } else if (isdigit(line.at(0))) {
            prime_5 = line.find('<') != std::string::npos;
            prime_3 = line.find('>') != std::string::npos;
            if (prime_5 && prime_3) {
                frame = FRAME_SELECTION_INTERNAL_FLAG;
            } else if (!prime_5 && !prime_3) {
                frame = FRAME_SELECTION_COMPLETE_FLAG;
            } else if (prime_5) {
                frame = FRAME_SELECTION_FIVE_FLAG;
            } else frame = FRAME_SELECTION_THREE_FLAG;
            pQuery_sequence = mpQueryData->get_sequence(seq_id);
            if (pQuery_sequence != nullptr) {
                pQuery_sequence->setFrame(frame);
            } else {
                throw ExceptionHandler("Sequence: " + seq_id + " not found in map during parsing of "
                                       "lst file at: " + mFinalLstPath,
                                       ERR_ENTAP_RUN_GENEMARK_PARSE);
            }
        }
    }
    FS_dprint("Success!");
}

/**
 * ======================================================================
 * Function void ModGeneMarkST::~ModGeneMarkST()
 *
 * Description          - Cleanup
 *
 * Notes                - Destructor
 *
 * @return              - None
 *
 * =====================================================================
 */
ModGeneMarkST::~ModGeneMarkST() {
    FS_dprint("Killing object - ModGeneMarkST");
}

/**
 * ======================================================================
 * Function void ModGeneMarkST::ModGeneMarkST()
 *
 * Description          - Initializes ModGeneMarkST object calling superclass AbstractFrame
 *                      - Sets final paths that GeneMarkS-T will create
 *
 * Notes                - Constructor
 *
 * @return              - None
 *
 * =====================================================================
 */
ModGeneMarkST::ModGeneMarkST(std::string &execution_stage_path, std::string &in_hits,
                             EntapDataPtrs &entap_data) :
    AbstractFrame(execution_stage_path, in_hits, entap_data, "GeneMarkS-T", DEFAULT_HEADERS) {
    mTranscriptomeFilename = mpFileSystem->get_filename(in_hits, true);
    mExePath = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_GENEMARKST_EXE);

    // Initialize FINAL output file paths
    // GenemarkS-T prints to CWD, they will be moved to these paths after execution
    mFinalFaaPath = PATHS(mModOutDir, mTranscriptomeFilename) + FileSystem::EXT_FAA;
    mFinalFnnPath = PATHS(mModOutDir, mTranscriptomeFilename) + FileSystem::EXT_FNN;
    mFinalLstPath = PATHS(mModOutDir, mTranscriptomeFilename + FileSystem::EXT_LST);
    mFinalGmstLogPath = PATHS(mModOutDir, GENEMARK_LOG_FILE);
    mFinalHmmPath = PATHS(mModOutDir, GENEMARK_HMM_FILE);
}

bool ModGeneMarkST::set_version() {
    return false;
}

bool ModGeneMarkST::is_executable(std::string &exe) {
    TerminalData terminalData;
    int return_code;

    terminalData.command = exe + " --version";
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;

    return_code = TC_execute_cmd(terminalData);  // GeneMarkS-T returns 1 with this command
    return_code >>= 8;                           // WARNING perl terminal command needs to be shifted by 8

    return return_code == GENEMARK_RETURN_OK;
}
