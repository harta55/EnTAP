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

#include "ModTransdecoder.h"
#include "../ExceptionHandler.h"

std::vector<ENTAP_HEADERS> ModTransdecoder::DEFAULT_HEADERS = {
        ENTAP_HEADER_FRAME
};

/**
 * ======================================================================
 * Function ModTransdecoder::ModTransdecoder(std::string &execution_stage_path, std::string &in_hits, EntapDataPtrs &entap_data,
                                 std::string &long_orfs_exe, std::string &predict_exe)
 *
 * Description          - Initialize ModTransdecoder object calling supers
 *                        (AbstractFrame and EntapModule)
 *
 * Notes                - Constructor
 *
 * @param execution_stage_path - Absolute path to Frame Selection directory
 * @param in_hits              - Input transcriptome that we would like to frame select
 * @param entap_data           - Pointers to necessary execution data
 * @param long_orfs_exe        - User input execution method for TransDecoder.LongOrfs
 * @param predict_exe          - User input execution method for TransDecoder.Predict
 *
 *
 * @return              - ModTransdecoder object
 *
 * =====================================================================
 */
ModTransdecoder::ModTransdecoder(std::string &execution_stage_path, std::string &in_hits, EntapDataPtrs &entap_data) :
    AbstractFrame(execution_stage_path, in_hits, entap_data, "TransDecoder", DEFAULT_HEADERS){
    std::string temp_filename;

    FS_dprint("Spawn Object - ModTransdecoder");

    mTransdecoderLongOrfsExe = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_TRANS_LONGORF_EXE);
    mTransdecoderPredictExe  = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_TRANS_PREDICT_EXE);
    mIsNoRefineStarts        = mpUserInput->has_input(INPUT_FLAG_TRANS_NO_REFINE_STARTS);
    mExePath = mTransdecoderLongOrfsExe;

    mMinProteinLength = mpUserInput->get_user_input<ent_input_uint_t >(INPUT_FLAG_TRANS_MIN_PROTEIN);

    // Determine TransDecoder output file naming
    temp_filename = mpFileSystem->get_filename(in_hits, true) + FILE_TRANSDECODER_SUFFIX + FileSystem::EXT_CDS;
    mOutputCDSFilePath = PATHS(mModOutDir, temp_filename);

    temp_filename = mpFileSystem->get_filename(in_hits, true) + FILE_TRANSDECODER_SUFFIX + FileSystem::EXT_PEP;
    mFinalFaaPath = PATHS(mModOutDir, temp_filename);   // PEP path
}

/**
 * ======================================================================
 * Function ModTransdecoder::~ModTransdecoder()
 *
 * Description          - Cleanup ModTransdecoder object
 *
 * Notes                - Destructor
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
ModTransdecoder::~ModTransdecoder() {
    FS_dprint("Killing object - ModTransdecoder");
}

/**
 * ======================================================================
 * Function EntapModule::ModVerifyData ModTransdecoder::verify_files()
 *
 * Description          - Checks whether we have previously ran TransDecoder
 *                      - Checks both PEP and CDS outputs to ensure we have them
 *
 * Notes                - output_paths of ModVerifyData is not currently referenced
 *
 *
 * @return              - ModVerifyData struct containing TRUE if we found
 *                        the files from a previous execution
 *
 * =====================================================================
 */
EntapModule::ModVerifyData ModTransdecoder::verify_files() {
    FS_dprint("Checking if TransDecoder files exist from a previous run...");

    ModVerifyData verify_data = ModVerifyData();

    verify_data.files_exist = true;

    // Only CDS and PEP files are required for continuation
    if (!mpFileSystem->file_exists(mOutputCDSFilePath) || !mpFileSystem->file_exists(mFinalFaaPath)) {
        FS_dprint("Missing either CDS or PEP files at:\nCDS: " + mOutputCDSFilePath + "\nPEP: " + mFinalFaaPath);
        verify_data.files_exist = false;
        // Remove both if exist
        mpFileSystem->delete_file(mOutputCDSFilePath);
        mpFileSystem->delete_file(mFinalFaaPath);
    } else {
        // Only push the final FAA file (unused currently)
        verify_data.output_paths.push_back(mFinalFaaPath);
    }
    return verify_data;
}

/**
 * ======================================================================
 * Function void ModTransdecoder::execute()
 *
 * Description          - Executes TransDecoder.LongOrfs software to train data
 *                      - Executes TransDecoder.Predict to predict ORFs
 *                      - Moves output from these (default is CWD) to
 *                        ModOutDir
 *
 * Notes                - WARNING TransDecoder prints to the CWD by default,
 *                        this will move some necessary files to the proper
 *                        directory, but not all. The "-O" flag added to later
 *                        versions does NOT work properly at this time
 *                      - Throws ExceptionHandler on errors
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModTransdecoder::execute() {
    std::string err_msg="";     // Message printed to user if TransDecoder fails
    std::string temp_filepath;
    std::string temp_filename;

    // Train data through Transdecoder.LongOrfs executable
    FS_dprint("Training TransDecoder data...");
    if (train_data(err_msg) != 0) {
        throw ExceptionHandler("Error in running TransDecoder.LongOrfs\nTransDecoder Error: " +
                err_msg, ERR_ENTAP_RUN_TRANSDECODER);
    }
    FS_dprint("Success!");

    // Predict coding regions through TransDecoder.Predict executable
    FS_dprint("Predicting coding regions through TransDecoder...");
    if (predict_frame(err_msg) != 0) {
        throw ExceptionHandler("Error in running TransDecoder.Predict\nTransDecoder Error: " +
                err_msg, ERR_ENTAP_RUN_TRANSDECODER);
    }

    FS_dprint("Success! Moving necessary output files...");
    try {
        // Move CDS file
        temp_filename = mpFileSystem->get_filename(mInputTranscriptome, true) + FILE_TRANSDECODER_SUFFIX + FileSystem::EXT_CDS;
        temp_filepath = PATHS(mpFileSystem->get_cur_dir(), temp_filename);
        mpFileSystem->rename_file(temp_filename, mOutputCDSFilePath);

        // Move PEP file
        temp_filename = mpFileSystem->get_filename(mInputTranscriptome, true) + FILE_TRANSDECODER_SUFFIX + FileSystem::EXT_PEP;
        temp_filepath = PATHS(mpFileSystem->get_cur_dir(), temp_filename);
        mpFileSystem->rename_file(temp_filename, mFinalFaaPath);

    } catch (const std::exception &err) {
        throw ExceptionHandler(err.what(), ERR_ENTAP_RUN_TRANSDECODER_MOVE);
    }
}

/**
 * ======================================================================
 * Function void ModTransdecoder::parse()
 *
 * Description          - Parses TransDecoder results, updating QueryData
 *                        structure and calculating statistics
 *
 * Notes                - None
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModTransdecoder::parse() {
    uint16 file_status;

    FS_dprint("Beginning to parse TransDecoder output...");

    // Ensure the files we need exist and are valid
    file_status = mpFileSystem->get_file_status(mFinalFaaPath);
    if (file_status != 0) {
        throw ExceptionHandler("TransDecoder PEP file invalid:\n" + mpFileSystem->print_file_status(file_status, mFinalFaaPath),
            ERR_ENTAP_RUN_TRANSDECODER_PARSE);
    }
    file_status = mpFileSystem->get_file_status(mOutputCDSFilePath);
    if (file_status != 0) {
        throw ExceptionHandler("TransDecoder CDS file invalid:\n" + mpFileSystem->print_file_status(file_status, mOutputCDSFilePath),
                               ERR_ENTAP_RUN_TRANSDECODER_PARSE);
    }
    // --------------- Files valid beyond this point ---------------

    try {
        // Pull relevant nucleotide, protein, and frame information and update QueryData
        parse_transdecoder_fasta(mFinalFaaPath, FileSystem::ENT_FILE_FASTA_FAA);
        parse_transdecoder_fasta(mOutputCDSFilePath, FileSystem::ENT_FILE_FASTA_FNN);

        FS_dprint("Success! TransDecoder files parsed");

        frame_calculate_statistics();

    } catch (const ExceptionHandler&e) {
        throw e;
    }

    FS_dprint("Success! TransDecoder data parsed and stats calculated");
}

/**
 * ======================================================================
 * Function bool ModDiamond::train_data(std::string &err_msg)
 *
 * Description          - Executes TransDecoder.LongOrfs software to train data
 *
 * Notes                - Runs to completion
 *                      - WARNING prints to the CWD
 *
 * @param err_msg       - Error message from TransDecoder execution (if exists)
 *
 *
 * @return              - INT status of TransDecoder.LongOrfs execution
 *
 * =====================================================================
 */
int ModTransdecoder::train_data(std::string &err_msg) {
    TerminalData  terminal_data;
    command_map_t tc_command_map;
    int ret=0;

    tc_command_map.emplace(CMD_TRANSCRIPTOME_INPUT, mInputTranscriptome);    // Transcriptome
    tc_command_map.emplace(CMD_MIN_PROTEIN_LENGTH, std::to_string(mMinProteinLength)); // Optional minimum protien length
//    tc_command_map.emplace(CMD_OUTPUT_DIR, mModOutDir); // Add output directory, DOESN't WORK!!! with latest version

    terminal_data.command = TC_generate_command(tc_command_map, mTransdecoderLongOrfsExe);
    terminal_data.base_std_path = PATHS(mModOutDir, STD_OUTPUT_TRAINING);
    terminal_data.print_files = true;
    terminal_data.suppress_std_err = true;  // Transdecoder goes ham on logging to err std, suppress it to not flood log

    ret = TC_execute_cmd(terminal_data);
    if (ret != 0) {
        err_msg = terminal_data.err_stream;
    }
    return ret;
}

/**
 * ======================================================================
 * Function bool ModDiamond::predict_frame(std::string &err_msg)
 *
 * Description          - Executes TransDecoder.Predict software to predict
 *                        protein reading frames
 *
 * Notes                - Runs to completion
 *                      - WARNING prints to the CWD
 *
 * @param err_msg       - Error message from TransDecoder execution (if exists)
 *
 *
 * @return              - INT status of TransDecoder.Predict execution
 *
 * =====================================================================
 */
int ModTransdecoder::predict_frame(std::string &err_msg) {
    TerminalData terminal_data;
    command_map_t tc_command_map;
    int ret;

    tc_command_map.emplace(CMD_TRANSCRIPTOME_INPUT, mInputTranscriptome);    // Transcriptome
    tc_command_map.emplace(CMD_SINGLE_BEST_ONLY, TC_NULL_ARGUMENT); // retain only one best ORF per sequence
//    tc_command_map.emplace(CMD_OUTPUT_DIR, mModOutDir); // Add output directory, DOESN't WORK!!! with latest version
    if (mIsNoRefineStarts) {
        tc_command_map.emplace(CMD_NO_REFINE_STARTS, TC_NULL_ARGUMENT);
    }

    terminal_data.command = TC_generate_command(tc_command_map, mTransdecoderPredictExe);
    terminal_data.base_std_path = PATHS(mModOutDir, STD_OUTPUT_PREDICTION);
    terminal_data.print_files = true;
    terminal_data.suppress_std_err = true; // Transdecoder goes ham on logging to err std, suppress it to not flood log

    ret = TC_execute_cmd(terminal_data);
    if (ret != 0) {
        err_msg = terminal_data.err_stream;
    }
    return ret;
}

/**
* ======================================================================
* Function bool ModTransdecoder::is_executable(std::string &long_orfs_exe,
 *                                              std::string &predict_exe)
*
* Description          - Determines whether both portions of TransDecoder
*                        (LongOrfs and Predict) will execute with the given
*                        execution paths
*
* Notes                - None
*
* @param long_orfs_exe - User input method of executing LongOrfs portion of TransDecoder
* @param predict_exe   - User input method of executing Predict portion of TransDecoder
*
*
* @return              - TRUE if TransDecoder is executable, FALSE otherwise
*
* =====================================================================
*/
bool ModTransdecoder::is_executable(std::string &long_orfs_exe, std::string &predict_exe) {
    TerminalData terminal_data;

    /* Test TransDecoder.LongOrfs */
    terminal_data.command = long_orfs_exe + " --version";
    terminal_data.print_files = false;
    terminal_data.suppress_std_err = false;

    FS_dprint("Testing TransDecoder.LongOrfs executable...");
    if (TC_execute_cmd(terminal_data) != 0) {
        FS_dprint("ERROR invalid TransDecoder.LongOrfs execution");
        return false;   // RETURN error in test
    }

    terminal_data = {};

    /* Test TransDecoder.Predict */
    terminal_data.command = predict_exe + " --version";
    terminal_data.print_files = false;

    FS_dprint("Testing TransDecoder.Predict executable...");
    if (TC_execute_cmd(terminal_data) != 0) {
        FS_dprint("ERROR invalid TransDecoder.Predict execution");
        return false;   // RETURN error in test
    }

    FS_dprint("Success! TransDecoder paths verified");
    return true;
}

/**
* ======================================================================
* Function void ModTransdecoder::parse_transdecoder_fasta(std::string &fasta_path,
 *                                          FileSystem::ENT_FILE_TYPES file_type)
*
* Description          - Parses TransDecoder CDS and PEP files and updates QueryData
*                        with Nucleotide/Protein/Frame information
*
* Notes                - Updates mpQueryData
*                      - Throws ExceptionHandler error on failure
*
* @param fasta_path    - Absolute path to fasta formatted file to parse
* @param file_type     - Type of FASTA file (FAA or FNN)
*
* @return              - None
*
* =====================================================================
*/
void ModTransdecoder::parse_transdecoder_fasta(std::string &fasta_path, FileSystem::ENT_FILE_TYPES file_type) {
    std::string                              line;
    std::string                              sequence;
    std::string                              seq_id;
    std::string                              frame;
    fp32                                     score;
    QuerySequence                            *pQuery_sequence=nullptr;

    std::ifstream in_file(fasta_path);
    try {
        while (true) {
            std::getline(in_file, line);
            if (line.empty() && !in_file.eof()) continue;
            // If line is the Sequence ID/Header line (signified by '>') or it is the last line of the file
            if (line.find(FileSystem::FASTA_FLAG) == 0 || in_file.eof()) {
                // Yes, this is a Sequence ID/Header line
                if (pQuery_sequence != nullptr) {
                    // If this is the last line of the file, add that sequence information
                    if (in_file.eof()) {
                        sequence += line + "\n";
                    }

                    // Finished collecting sequence information, update QueryData
                    // In the event that there are multiple frames, ensure we have the best one
                    if (score >= pQuery_sequence->getMFrameScore()) {
                        if (file_type == FileSystem::ENT_FILE_FASTA_FAA) {
                            pQuery_sequence->set_sequence_p(sequence);
                        } else if (file_type == FileSystem::ENT_FILE_FASTA_FNN) {
                            pQuery_sequence->set_sequence_n(sequence);
                        }
                        pQuery_sequence->setMFrameScore(score);
                        pQuery_sequence->setFrame(frame);
                    }
                }

                if (in_file.eof()) {
                    // EXIT WHILE if we have reached the end of the file
                    break;
                } else {
                    // Set the beginning of this new sequence to the Sequence ID/Header
                    parse_transdecoder_fasta_header(seq_id, line, frame, score);
                    pQuery_sequence = mpQueryData->get_sequence(seq_id);
                    if (pQuery_sequence == nullptr) {
                        throw ExceptionHandler("Unable to find correct sequence ID from TransDecoder line: " + line,
                            ERR_ENTAP_RUN_TRANSDECODER_PARSE);
                    } else {
                        // Ensure we have the correct ID
                        seq_id = pQuery_sequence->getMSequenceID();
                        sequence = FileSystem::FASTA_FLAG + seq_id + "\n";

                        // Verify frame is correct and update with standard
                        verify_ORF(frame);
                    }
                }

            } else {
                // No, must be a generic sequence line

                // TransDecoder may add '*' at the end of a sequence. Delete this in order to prevent
                // incompatibility with other SW
                if (line.back() =='*') {
                    line.pop_back();
                }
                sequence += line + "\n";
            }
        }
        in_file.close();
    } catch (const std::exception &e) {
        throw ExceptionHandler("Error in TransDecoder parsing of file at: " + fasta_path + "\nError:\n" +
                e.what(), ERR_ENTAP_RUN_TRANSDECODER_PARSE);
    }
}

/**
* ======================================================================
* Function bool ModTransdecoder::verify_ORF(std::string &frame)
*
* Description          - Ensures the ORF pulled from the TransDecoder header
*                        is valid and can be mapped to a known ORF
*
* Notes                - Sets frame to mapped frame
*
* @param frame    - Set to ORF string mapped from TransDecoder string
*
* @return              - None
*
* =====================================================================
*/
bool ModTransdecoder::verify_ORF(std::string &frame) {
    bool ret = true;

    if (frame == TRANSDECODER_FLAG_COMPLETE_FRAME) {
        frame = FRAME_SELECTION_COMPLETE_FLAG;
    } else if (frame == TRANSDECODER_FLAG_INTERNAL_FRAME) {
        frame = FRAME_SELECTION_INTERNAL_FLAG;
    } else if (frame == TRANSDECODER_FLAG_3PRIME_FRAME) {
        frame = FRAME_SELECTION_THREE_FLAG;
    } else if (frame == TRANSDECODER_FLAG_5PRIME_FRAME) {
        frame = FRAME_SELECTION_FIVE_FLAG;
    } else {
        throw ExceptionHandler("Frame could not be pulled from TransDecoder output: " + frame,
            ERR_ENTAP_RUN_TRANSDECODER_PARSE);
    }
    return ret;
}

/**
* ======================================================================
* Function std::string ModTransdecoder::format_fasta_header(std::string &seq_id,
 *                                                          std::string &line, std::string &frame)
*
* Description          - This routine will parse the TransDecoder fasta
*                        header (CDS/PEP files) and set the sequence ID (seq_id)
*                        and frame
*
*
* Notes                - Sets seq_id to the sequence ID
*                      - No std C++11 REGEX for older compiler support
*
* @param seq_id        - This value is set to the sequence ID from the parsed header
* @param line          - Entire header line (starting with '>')
* @param frame         - This value is set to the ORF type pulled from the header
* @param score         - This value is set ot the score parsed from the header
*
* @return              - Sets seq_id and frame
*
* =====================================================================
*/
void ModTransdecoder::parse_transdecoder_fasta_header(std::string &seq_id, std::string &line, std::string &frame, fp32 &score) {
    uint64 ind1;
    uint64 ind2;
    std::string ret;

    // Find sequence ID from entire header
    // Typically, this is the format we are looking for in one line:
    //      >TR.ECb10001|c0_g1_i1.p1 GENE.TR.ECb10001|c0_g1_i1~~TR.ECb10001|c0_g1_i1.p1
    //      ORF type:internal len:110 (+),score=0.03 TR.ECb10001|c0_g1_i1:1-327(+)
    ind1 = line.find(FLAG_QUERYID_BEGIN);
    ind2 = line.find(FLAG_QUERYID_END);
    if (ind1 == std::string::npos || ind2 == std::string::npos) {
        // WARNING there is some interaction between Transdecoder and certain transcriptomes that report the following
        //  format in the .pep file, this must be parsed differently:
        //      >TRINITY_DN0_c1_g1_i3.p1 TRINITY_DN0_c1_g1~~TRINITY_DN0_c1_g1_i3.p1
        //      ORF type:internal len:121 (-),score=46.75 TRINITY_DN0_c1_g1_i3:3-362(-)
        uint64 temp = line.find(FLAG_SCORE_BEGIN);
        ind1 = line.find(FLAG_QUERYID_V2_BEGIN, temp);
        ind2 = line.find(FLAG_QUERYID_V2_END, ind1);
        if (ind1 == std::string::npos || ind2 == std::string::npos) {
            // If even this doesnt work, throw error
            FS_dprint("Unable to find sequence id from header: " + line);
            throw ExceptionHandler("Unable to find sequence id from header: " + line, ERR_ENTAP_RUN_TRANSDECODER_PARSE);
        } else {
            seq_id = line.substr(ind1 + 1, ind2 - ind1 - 1);
        }

    } else {
        seq_id = line.substr(ind1 + FLAG_QUERYID_BEGIN.size(), ind2 - ind1 - FLAG_QUERYID_BEGIN.size());
    }

    // Find frame from entire header
    ind1 = line.find(FLAG_FRAME_BEGIN);
    ind2 = line.find(FLAG_FRAME_END, ind1);
    if (ind1 == std::string::npos || ind2 == std::string::npos) {
        throw ExceptionHandler("Unable to find sequence ORF from header: " + line, ERR_ENTAP_RUN_TRANSDECODER_PARSE);
    } else {
        frame = line.substr(ind1 + FLAG_FRAME_BEGIN.size(), ind2 - ind1 - FLAG_FRAME_BEGIN.size());
    }

    // Find score from entire header
    ind1 = line.find(FLAG_SCORE_BEGIN);
    ind2 = line.find(FLAG_SCORE_END, ind1);
    if (ind1 == std::string::npos || ind2 == std::string::npos) {
        // Do not error out for this
        score = 0.0;
    } else {
        std::string score_str = line.substr(ind1 + FLAG_SCORE_BEGIN.size(), ind2 - ind1 - FLAG_SCORE_BEGIN.size());
        score = std::stof(score_str);
    }
}

bool ModTransdecoder::set_version() {
    return false;
}
