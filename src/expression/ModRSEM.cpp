/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
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

//**************************************************************


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
std::pair<bool, std::string> ModRSEM::verify_files() {

    boostFS::path                   file_name(_inpath);

    while (file_name.has_extension()) file_name = file_name.stem();
    _filename = file_name.string();
    _exp_out = PATHS(_expression_outpath, _filename);
    _rsem_out = _exp_out + RSEM_OUT_FILE;
    if (FS_file_exists(_rsem_out)) {
        FS_dprint("File found at " + _rsem_out +  "\nmoving to filter transcriptome");
        return std::make_pair(true, "");
    }
    FS_dprint("File not found at " + _rsem_out +  " Continuing RSEM run.");

    return std::make_pair(false, "");
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

    boostFS::path                   bam_ext(_alignpath);
    std::string                     bam;
    std::string                     rsem_arg;
    std::string                     out_path;
    std::string                     ref_path;   // Reference path for RSEM
    std::string                     ref_exe;
    std::string                     expression_exe;
    std::string                     std_out;    // Outpath for std err and out

    bam_ext = bam_ext.extension();
    bam = bam_ext.string();
    std::transform(bam.begin(), bam.end(), bam.begin(), ::tolower);
    bam = bam.substr(1,3); // just extract extension
    //todo separate into methods
    if (!rsem_validate_file(_filename)){
        throw ExceptionHandler("Alignment file can't be validated by RSEM. Check error file!",
                               ENTAP_ERR::E_RUN_RSEM_VALIDATE);
    }
    // Now have valid BAM file to run rsem
    FS_dprint("Alignment file valid. Preparing reference...");

    ref_exe  = PATHS(_exe_path, RSEM_PREP_REF_EXE);
    ref_path = PATHS(_expression_outpath, _filename) + "_ref";

    // Prepare reference command
    rsem_arg = ref_exe + " "
               + _inpath + " "
               + ref_path;
    std_out = PATHS(_expression_outpath, _filename) + STD_REF_OUT;
    FS_dprint("Executing following command\n" + rsem_arg);
    execute_cmd(rsem_arg.c_str(), std_out);
    FS_dprint("Reference successfully created");

    FS_dprint("Running expression analysis...");
    expression_exe = PATHS(_exe_path, RSEM_CALC_EXP_EXE);
    rsem_arg = expression_exe +
               " --" + bam +
               " -p " + std::to_string(_threads) + " " +
               _alignpath +" "+
               ref_path + " " +
               _exp_out;
    if (!_issingle) rsem_arg += " --paired-end";
    std_out = PATHS(_expression_outpath, _filename) + STD_EXP_OUT;
    FS_dprint("Executing following command\n" + rsem_arg);
    if (execute_cmd(rsem_arg.c_str(), std_out)!=0) {
        throw ExceptionHandler("Error in running expression analysis",ENTAP_ERR::E_INIT_TAX_READ);
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
std::string ModRSEM::filter() {
    FS_dprint("Beginning to filter transcriptome...");

    uint32              count_removed=0;
    uint32              count_kept=0;
    uint32              count_total=0;      // Used to warn user if high percentage is removed
    uint32              min_removed=10000;
    uint32              min_selected=10000;
    uint32              max_removed=10000;
    uint32              max_selected=10000;
    uint64              total_removed_len=0;
    uint64              total_kept_len=0;
    uint16              length;
    fp32                in_len;
    fp32                e_leng;
    fp32                e_count;
    fp32                tpm;
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
    std::string         fig_txt_box_path;
    std::string         fig_png_box_path;
    std::string         min_kept_seq;
    std::string         max_kept_seq;
    std::string         max_removed_seq;
    std::string         min_removed_seq;
    std::stringstream   out_msg;
    std::vector<uint16> all_kept_lengths;
    std::vector<uint16> all_lost_lengths;
    std::pair<uint64,uint64> kept_n;
    GraphingData        graphingStruct;
    QUERY_MAP_T         *MAP;

    MAP = pQUERY_DATA->get_sequences_ptr();

    if (!FS_file_exists(_rsem_out)) {
        throw ExceptionHandler("File does not exist at: " + _rsem_out, ENTAP_ERR::E_RUN_RSEM_EXPRESSION);
    }

    FS_delete_dir(_processed_path);
    FS_create_dir(_processed_path);
    FS_create_dir(_figure_path);
    boostFS::path path (_rsem_out);

    fig_txt_box_path = PATHS(_figure_path, GRAPH_TXT_BOX_PLOT);
    fig_png_box_path = PATHS(_figure_path, GRAPH_PNG_BOX_PLOT);

    std::ofstream file_fig_box(fig_txt_box_path, std::ios::out | std::ios::app);
    file_fig_box << "flag\tsequence length" << std::endl;    // First line placeholder, not used

    while (path.has_extension()) path = path = path.stem();
    original_filename = path.string();
    removed_filename  = original_filename + RSEM_OUT_REMOVED;
    kept_filename     = original_filename + RSEM_OUT_KEPT;
    out_kept    = PATHS(_processed_path, kept_filename);
    out_removed = PATHS(_processed_path, removed_filename);
    std::ofstream out_file(out_kept, std::ios::out | std::ios::app);
    std::ofstream removed_file(out_removed, std::ios::out | std::ios::app);
    io::CSVReader<RSEM_COL_NUM, io::trim_chars<' '>,
    io::no_quote_escape<'\t'>> in(_rsem_out);
    in.next_line();
    while (in.read_row(geneid, transid, in_len, e_leng, e_count, tpm, fpkm_val)) {
        count_total++;
        pQUERY_DATA->trim_sequence_header(geneid,geneid);
        QUERY_MAP_T::iterator it = MAP->find(geneid);
        if (it == MAP->end()) {
            throw ExceptionHandler("Unable to find sequence: " + geneid,
                                   ENTAP_ERR::E_RUN_RSEM_EXPRESSION_PARSE);
        }
        QuerySequence *querySequence = it->second;
        length = (uint16)querySequence->getSeq_length();
        if (fpkm_val > _fpkm) {
            // Kept sequence
            out_file << querySequence->get_sequence() << std::endl;
            querySequence->set_fpkm(fpkm_val);
            file_fig_box << GRAPH_KEPT_FLAG << '\t' << std::to_string(length) << std::endl;
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
            file_fig_box << GRAPH_REJECTED_FLAG << '\t' << std::to_string(length) << std::endl;

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
    out_msg<<std::fixed<<std::setprecision(2);
    out_msg <<
            ENTAP_STATS::SOFTWARE_BREAK     <<
            "Expression Filtering - RSEM\n" <<
            ENTAP_STATS::SOFTWARE_BREAK     <<
            "Total sequences kept: "        << count_kept     <<
            "\nTotal sequences removed: "   << count_removed;


    if (count_kept > 0) {
        rejected_percent = ((fp32)count_removed / count_kept) * 100;
        avg_kept = (fp32) total_kept_len / count_kept;
        kept_n = pQUERY_DATA->calculate_N_vals(all_kept_lengths, total_kept_len);
        out_msg <<
                ENTAP_STATS::SOFTWARE_BREAK                                     <<
                "Expression Filtering: New Reference Transcriptome Statistics"  <<
                ENTAP_STATS::SOFTWARE_BREAK                                     <<
                "\nTotal sequenes: "                    << count_kept     <<
                "\nTotal length of transcriptome (bp)"  << total_kept_len <<
                "\nAverage length (bp): "               << avg_kept       <<
                "\nn50: "                               << kept_n.first   <<
                "\nn90: "                               << kept_n.second  <<
                "\nLongest sequence (bp): " << max_selected << " (" << max_kept_seq << ")" <<
                "\nShortest sequence (bp): "<< min_selected << " (" << min_kept_seq << ")\n";
    } else {
        throw ExceptionHandler("Error in filtering transcriptome, no sequences kept",
                               ENTAP_ERR::E_RUN_RSEM_EXPRESSION);
    }

    if (count_removed > 0) {
        avg_removed = (fp32) total_removed_len / count_removed;
        std::pair<uint64, uint64> removed_n =
                pQUERY_DATA->calculate_N_vals(all_lost_lengths,total_removed_len);
        out_msg <<
                "\nRemoved Sequences (no frame):"       <<
                "\nTotal sequences: "                     << count_removed    <<
                "\nAverage sequence length(bp): "         << avg_removed      <<
                "\nn50: "                                 << removed_n.first  <<
                "\nn90: "                                 << removed_n.second <<
                "\nLongest sequence(bp): "  << max_removed<< " (" << max_removed_seq << ")" <<
                "\nShortest sequence(bp): " << min_removed<< " (" << min_removed_seq << ")" <<"\n";

        if (rejected_percent > REJECTED_ERROR_CUTOFF) {
            // Warn user high percentage of transcriptome was rejected
            out_msg << "\nWarning: A high percentage of the transcriptome was removed: " << rejected_percent;
        }
    }

    out_str = out_msg.str();
    FS_print_stats(out_str);
    out_file.close();
    removed_file.close();
    file_fig_box.close();
    FS_dprint("Success!");
    //--------------------------------------------------------//


    //------------------------Graphing------------------------//
    FS_dprint("Beginning to send data to graphing manager...");
    graphingStruct.text_file_path   = fig_txt_box_path;
    graphingStruct.graph_title      = GRAPH_TITLE_BOX_PLOT;
    graphingStruct.fig_out_path     = fig_png_box_path;
    graphingStruct.software_flag    = GRAPH_EXPRESSION_FLAG;
    graphingStruct.graph_type       = GRAPH_BOX_FLAG;
    pGraphingManager->graph(graphingStruct);
    FS_dprint("Success!");
    //--------------------------------------------------------//

    return out_kept;
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

    rsem_arg = PATHS(_exe_path, RSEM_SAM_VALID) + " " + _alignpath;
    out_path = PATHS(_expression_outpath, filename) + STD_VALID_OUT;
    // only thrown in failure in calling rsem
    FS_dprint("Executing RSEM command:\n" + rsem_arg);
    if (execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) return false;
    FS_dprint("RSEM validate executed successfully");
    // RSEM does not always return error code if file is invalid, only seen in .err
    return (FS_file_no_lines(out_path + EXT_ERR));
}


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

    bam_out  = PATHS(_expression_outpath, file_name);
    rsem_arg = PATHS(_exe_path, RSEM_CONV_SAM) +
               " -p " + std::to_string(_threads) + " " + _alignpath + " " + bam_out;
    std_out  = _expression_outpath + file_name + STD_CONVERT_SAM;
    if (execute_cmd(rsem_arg.c_str(), std_out)!=0)return false;
    if (!FS_file_no_lines(std_out+EXT_ERR)) return false;
    _alignpath = bam_out + EXT_BAM;
    return true;
}


/**
 * ======================================================================
 * Function void ModRSEM::set_data(int thread, float fpmk, bool paired)
 *
 * Description          - Sets specific data used by each module (probably
 *                        will change in adding modules)
 *
 * Notes                - None
 *
 * @param thread        - Thead count
 * @param fpkm          - User selected FPKM threshold (for filtering)
 * @param paired        - Yes/no if reads were paired-end during alignment
 *
 * @return              - True if file is empty
 *
 * =====================================================================
 */
void ModRSEM::set_data(int thread, float fpkm, bool single) {
    _threads = thread;
    _fpkm = fpkm;
    _issingle = single;
}


