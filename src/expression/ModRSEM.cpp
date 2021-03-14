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


//*********************** Includes *****************************
#include "ModRSEM.h"
#include "../TerminalCommands.h"
#include "../QueryData.h"

//**************************************************************

std::vector<ENTAP_HEADERS> ModRSEM::DEFAULT_HEADERS = {
        ENTAP_HEADER_EXP_FPKM,
        ENTAP_HEADER_EXP_TPM,
        ENTAP_HEADER_EXP_E_LENGTH
};


ModRSEM::ModRSEM(std::string &execution_stage_path, std::string &in_hits, EntapDataPtrs &entap_data,
                 std::string &align) :
        AbstractExpression(execution_stage_path, in_hits, entap_data, "RSEM", DEFAULT_HEADERS){
    FS_dprint("Spawn Object - ModRSEM");

    mCalcExpressionExe = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_RSEM_CALC_EXPRES);
    mSamValidExe       = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_RSEM_SAM_VALID);
    mPrepReferenceExe  = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_RSEM_PREP_REF);
    mConvertSamExe     = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_RSEM_CONVERT_SAM);
    mExePath           = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_RSEM_CALC_EXPRES);
}

ModRSEM::~ModRSEM() {
    FS_dprint("Killing object - ModRSEM");
}


/**
 * ======================================================================
 * Function std::pair<bool, std::string> ModRSEM::verify_files()
 *
 * Description          - Verifies whether the user has ran RSEM previously
 *                      - Checks for specific output file
 *
 * Notes                - None
 *
 * @return              - Pair:
 *                              - True if file has been found
 *                              - Path to the file (not used, carried over)
 *
 * =====================================================================
 */
EntapModule::ModVerifyData ModRSEM::verify_files() {
    ModVerifyData modVerifyData;

    mFilename = mpFileSystem->get_filename(mInputTranscriptome, false);
    mExpressionOut  = PATHS(mModOutDir, mFilename);
    mRsemOut = mExpressionOut + RSEM_OUT_FILE;
    if (mpFileSystem->file_exists(mRsemOut)) {
        FS_dprint("File found at " + mRsemOut +  "\nmoving to filter transcriptome");
        modVerifyData.files_exist = true;
    } else {
        FS_dprint("File not found at " + mRsemOut +  " Continuing RSEM run.");
        modVerifyData.files_exist = false;
    }
    return modVerifyData;
}


/**
 * ======================================================================
 * Function ModRSEM::execute(std::map<std::string, QuerySequence> &)
 *
 * Description          - Responsible for executing RSEM with user
 *                        alignment file
 *                      - Verifies input SAM/BAM file through RSEM
 *                      - Prepares a reference from the transcriptome
 *                      - Execution expression analysis
 *
 * Notes                - Incorporates --paired-end flag
 *
 * @param map           - Not used
 *                        complete gene
 * @return              - None
 *
 * =====================================================================
 */
void ModRSEM::execute() {
    // return path
    FS_dprint("Running RSEM...");

    std::string                     bam;
    std::string                     rsem_arg;
    std::string                     out_path;
    std::string                     ref_path;   // Reference path for RSEM
    std::string                     ref_exe;
    std::string                     expression_exe;
    std::string                     std_out;    // Outpath for std err and out

    bam = mpFileSystem->get_file_extension(mAlignPath, true);
    LOWERCASE(bam);
    if (!rsem_validate_file(mFilename)){
        throw ExceptionHandler("Alignment file can't be validated by RSEM. Check error file!",
                               ERR_ENTAP_RUN_RSEM_VALIDATE);
    }
    // Now have valid BAM file to run rsem
    FS_dprint("Alignment file valid. Preparing reference...");

    // Prepare reference command
    if (rsem_generate_reference(ref_path)) {
        FS_dprint("Reference successfully created");
    } else {
        throw ExceptionHandler("Unable to generate RSEM reference", ERR_ENTAP_RUN_RSEM_EXPRESSION);
    }

    // Run expression analysis
    FS_dprint("Running expression analysis...");
    if (rsem_expression_analysis(ref_path, bam)) {
        FS_dprint("Expression analysis complete!");
    } else {
        throw ExceptionHandler("Error in running expression analysis",ERR_ENTAP_RUN_RSEM_EXPRESSION);
    }
}


/**
 * ======================================================================
 * Function std::string ModRSEM::filter(std::map<std::string, QuerySequence> & MAP)
 *
 * Description          - Handles filtering of transcriptome based on
 *                        user selected FPKM threshold
 *                      - Updates master query sequence map
 *
 * Notes                - None
 *
 * @param MAP           - Master QuerySequence map
 *
 * @return              - Output path of sequences that were kept
 *
 * =====================================================================
 */
void ModRSEM::parse() {
    FS_dprint("Beginning to filter transcriptome...");

    uint32              count_removed=0;
    uint32              count_kept=0;
    uint32              count_total=0;      // Used to warn user if high percentage is removed
    uint32              min_removed=0xFFFFFFFF;
    uint32              min_selected=0xFFFFFFFF;
    uint32              max_removed=0;
    uint32              max_selected=0;
    uint64              total_removed_len=0;
    uint64              total_kept_len=0;
    uint16              length;
    fp32                in_len; // Length
    fp32                e_leng; // Effective length
    fp32                e_count;
    fp64                tpm;
    fp32                fpkm_val;
    fp32                rejected_percent=0;
    fp32                avg_removed;
    fp32                avg_kept;
    std::string         geneid;
    std::string         transid;
    std::string         out_str;
    std::string         out_kept;
    std::string         out_removed;
    std::string         kept_filename;
    std::string         original_filename;
    std::string         removed_filename;
    std::string         min_kept_seq;
    std::string         max_kept_seq;
    std::string         max_removed_seq;
    std::string         min_removed_seq;
    std::stringstream   out_msg;
    std::vector<uint16> all_kept_lengths;
    std::vector<uint16> all_lost_lengths;
    std::pair<uint64,uint64> kept_n;

    if (!mpFileSystem->file_exists(mRsemOut)) {
        throw ExceptionHandler("File does not exist at: " + mRsemOut,
                               ERR_ENTAP_RUN_RSEM_EXPRESSION);
    }

    // Setup figure files, directories already created
    GraphingManager::GraphingData graph_box_plot;
    graph_box_plot.x_axis_label = "Flag";
    graph_box_plot.y_axis_label = "Sequence Length";
    graph_box_plot.text_file_path = PATHS(mFigureDir, GRAPH_TXT_BOX_PLOT);
    graph_box_plot.fig_out_path   = PATHS(mFigureDir, GRAPH_PNG_BOX_PLOT);
    graph_box_plot.graph_title    = GRAPH_TITLE_BOX_PLOT;
    graph_box_plot.graph_type     = GraphingManager::ENT_GRAPH_BOX_PLOT_VERTICAL;
    mpGraphingManager->initialize_graph_data(graph_box_plot);

    // Setup processed file paths, directories already created
    original_filename = mFilename;
    removed_filename  = original_filename + RSEM_OUT_REMOVED;
    kept_filename     = original_filename + RSEM_OUT_KEPT;
    out_kept    = PATHS(mProcDir, kept_filename);
    out_removed = PATHS(mProcDir, removed_filename);
    std::ofstream out_file(out_kept, std::ios::out | std::ios::app);
    std::ofstream removed_file(out_removed, std::ios::out | std::ios::app);

    // Begin to iterate through RSEM output file
    io::CSVReader<RSEM_COL_NUM, io::trim_chars<' '>,
    io::no_quote_escape<'\t'>> in(mRsemOut);
    in.next_line();
    while (in.read_row(geneid, transid, in_len, e_leng, e_count, tpm, fpkm_val)) {
        count_total++;
        mpQueryData->trim_sequence_header(geneid,geneid);
        QuerySequence *querySequence = mpQueryData->get_sequence(geneid);
        if (querySequence == nullptr) {
            throw ExceptionHandler("Unable to find sequence: " + geneid + " there may be a discrepancy between"
                                                                          " sequence headers in your transcriptome and "
                                                                          "headers in your BAM/SAM file. Try trimming"
                                                                          " your sequence headers to the first space and re-running.",
                                   ERR_ENTAP_RUN_RSEM_EXPRESSION_PARSE);
        }
        querySequence->set_fpkm(fpkm_val);
        querySequence->setMTPM(tpm);
        querySequence->setMEffectiveLength(e_leng);
        length = (uint16) querySequence->get_sequence_length();
        if (fpkm_val > mFPKM) {
            // Kept sequence
            out_file << querySequence->get_sequence() << std::endl;
            mpGraphingManager->add_datapoint(graph_box_plot.text_file_path, {GRAPH_KEPT_FLAG, std::to_string(length)});
            //TODO move to QueryData
            if (length < min_selected) {
                min_selected = length;
                min_kept_seq = geneid;
            }
            if (length > max_selected) {
                max_selected = length;
                max_kept_seq = geneid;
            }
            all_kept_lengths.push_back(length);
            total_kept_len += length;
            count_kept++;
        } else {
            // Removed sequence
            querySequence->QUERY_FLAG_CLEAR(QuerySequence::QUERY_EXPRESSION_KEPT);
            removed_file << querySequence->get_sequence() << std::endl;
            mpGraphingManager->add_datapoint(graph_box_plot.text_file_path, {GRAPH_REJECTED_FLAG, std::to_string(length)});

            if (length < min_removed) {
                min_removed = length;
                min_removed_seq = geneid;
            }
            if (length > max_removed) {
                max_removed_seq = geneid;
                max_removed = length;
            }
            all_lost_lengths.push_back(length);
            total_removed_len += length;
            count_removed++;
        }
    }
    FS_dprint("File successfully filtered. Outputs at:\n" + out_kept + " and:\n" + out_removed);

    //-----------------------STATISTICS-----------------------//
    FS_dprint("Beginning to calculate statistics...");
    mpFileSystem->format_stat_stream(out_msg, "Expression Filtering (RSEM) with FPKM Cutoff " + float_to_string(mFPKM));
    out_msg <<
            "Total sequences kept: "        << count_kept     <<
            "\nTotal sequences removed: "   << count_removed  <<std::endl;


    if (count_kept > 0) {
        rejected_percent = ((fp32)count_removed / count_total) * 100;
        avg_kept = (fp32) total_kept_len / count_kept;
        kept_n = mpQueryData->calculate_N_vals(all_kept_lengths, total_kept_len);
        mpFileSystem->format_stat_stream(out_msg, "Expression Filtering: New Reference Transcriptome Statistics");
        out_msg <<
                "\nTotal sequenes: "                    << count_kept     <<
                "\nTotal length of transcriptome (bp)"  << total_kept_len <<
                "\nAverage length (bp): "               << avg_kept       <<
                "\nn50: "                               << kept_n.first   <<
                "\nn90: "                               << kept_n.second  <<
                "\nLongest sequence (bp): " << max_selected << " (" << max_kept_seq << ")" <<
                "\nShortest sequence (bp): "<< min_selected << " (" << min_kept_seq << ")\n";
    } else {
        throw ExceptionHandler("Error in filtering transcriptome, no sequences kept",
                               ERR_ENTAP_RUN_RSEM_EXPRESSION);
    }

    if (count_removed > 0) {
        avg_removed = (fp32) total_removed_len / count_removed;
        std::pair<uint64, uint64> removed_n =
                mpQueryData->calculate_N_vals(all_lost_lengths,total_removed_len);
        out_msg <<
                "\nRemoved Sequences (under FPKM threshold):"       <<
                "\nTotal sequences: "                     << count_removed    <<
                "\nAverage sequence length(bp): "         << avg_removed      <<
                "\nn50: "                                 << removed_n.first  <<
                "\nn90: "                                 << removed_n.second <<
                "\nLongest sequence(bp): "  << max_removed<< " (" << max_removed_seq << ")" <<
                "\nShortest sequence(bp): " << min_removed<< " (" << min_removed_seq << ")" <<"\n";

        if (rejected_percent > REJECTED_ERROR_CUTOFF) {
            // Warn user high percentage of transcriptome was rejected
            out_msg << "\nWARNING: A high percentage of the transcriptome was removed: " << rejected_percent;
        }
    } else {
        out_msg << "\nWARNING: No sequences were removed from Expression Filtering";
    }

    out_str = out_msg.str();
    mpFileSystem->print_stats(out_str);
    out_file.close();
    removed_file.close();
    FS_dprint("Success!");
    //--------------------------------------------------------//


    //------------------------Graphing------------------------//
    FS_dprint("Beginning to send data to graphing manager...");
    mpGraphingManager->graph_data(graph_box_plot.text_file_path);
    FS_dprint("Success!");
    //--------------------------------------------------------//

    mFinalFasta = out_kept;
}


/**
 * ======================================================================
 * Function ModRSEM::rsem_validate_file(std::string filename)
 *
 * Description          - Executes rsem-sam-validator
 *                      - Determines if user inpuuted BAM/SAM file is valid
 *                        and will continue to expression analysis
 *
 * Notes                - None
 *
 * @param filename      - Transcriptome name just to label output
 *
 * @return              - True if a valid SAM/BAM file
 *
 * =====================================================================
 */
bool ModRSEM::rsem_validate_file(std::string filename) {
    FS_dprint("File is detected to be sam file, running validation");

    std::string rsem_arg;
    std::string out_path;
    TerminalData terminalData;

    rsem_arg = mSamValidExe + " " + mAlignPath;
    out_path = PATHS(mModOutDir, filename) + STD_VALID_OUT;

    terminalData.command        = rsem_arg;
    terminalData.print_files    = true;
    terminalData.suppress_std_err = false;
    terminalData.base_std_path  = PATHS(mModOutDir, filename) + STD_VALID_OUT;


    if (TC_execute_cmd(terminalData)!=0) {
        throw ExceptionHandler("Error in validating SAM file\nRSEM Error:\n" +
                               terminalData.err_stream, ERR_ENTAP_RUN_RSEM_EXPRESSION);
    }
    FS_dprint("RSEM validate executed successfully");
    // RSEM does not always return error code if file is invalid, only seen in .err
    return (mpFileSystem->file_no_lines(out_path + FileSystem::EXT_ERR));
}

#if 0
/**
 * ======================================================================
 * Function ModRSEM::rsem_conv_to_bam(std::string file_name)
 *
 * Description          - Executes convert-sam-for-rsem
 *                      - Converts user inputed file to BAM format
 *
 * Notes                - Not used, done automatically if needed
 *
 * @param filename      - Transcriptome name just to label output
 *
 * @return              - True if conversion was successful
 *
 * =====================================================================
 */
bool ModRSEM::rsem_conv_to_bam(std::string file_name) {
    FS_dprint("Converting SAM to BAM");

    std::string bam_out;
    std::string rsem_arg;
    std::string std_out;

    bam_out  = PATHS(mModOutDir, file_name);
    rsem_arg = PATHS(mExePath, RSEM_CONV_SAM) +
               " -p " + std::to_string(mThreads) + " " + mAlignPath + " " + bam_out;
    std_out  = mModOutDir + file_name + STD_CONVERT_SAM;
    if (TC_execute_cmd(rsem_arg.c_str(), std_out)!=0)return false;
    if (!pFileSystem->file_no_lines(std_out+FileSystem::EXT_ERR)) return false;
    mAlignPath = bam_out + FileSystem::EXT_BAM;
    return true;
}
#endif

bool ModRSEM::rsem_generate_reference(std::string& reference_path_out) {
    std::string       ref_path;   // reference path
    std::string       rsem_arg;
    int32             err_code;
    TerminalData      terminalData;

    ref_path = PATHS(mModOutDir, mFilename) + "_ref";
    rsem_arg = mPrepReferenceExe + " "
               + mInputTranscriptome + " "
               + ref_path;

    terminalData.command       = rsem_arg;
    terminalData.base_std_path = PATHS(mModOutDir, mFilename) + STD_REF_OUT;
    terminalData.print_files   = true;
    terminalData.suppress_std_err = false;

    reference_path_out = ref_path;

    // Execute command
    err_code = TC_execute_cmd(terminalData);

    if (err_code != 0) {
        throw ExceptionHandler("Error generating RSEM reference\nRSEM Error:\n" +
                               terminalData.err_stream, ERR_ENTAP_RUN_RSEM_EXPRESSION);
    }

    return err_code == 0;
}

bool ModRSEM::rsem_expression_analysis(std::string& ref_path, std::string& bam) {
    std::string rsem_arg;
    int32             err_code;
    TerminalData      terminalData;

    rsem_arg = mCalcExpressionExe +
               " --" + bam +
               " -p " + std::to_string(mThreads) + " " +
               mAlignPath +" "+
               ref_path + " " +
               mExpressionOut;
    if (!mIsSingle) rsem_arg += " --paired-end";

    terminalData.command        = rsem_arg;
    terminalData.print_files    = true;
    terminalData.suppress_std_err = false;
    terminalData.base_std_path  = PATHS(mModOutDir, mFilename) + STD_EXP_OUT;


    err_code = TC_execute_cmd(terminalData);

    if (err_code != 0) {
        throw ExceptionHandler("Error in running RSEM Expression Analysis\nRSEM Error:\n" +
                terminalData.err_stream, ERR_ENTAP_RUN_RSEM_EXPRESSION);
    }
    return err_code == 0;
}


std::string ModRSEM::get_final_fasta() {
    return this->mFinalFasta;
}

bool ModRSEM::set_version() {
    return false;
}
