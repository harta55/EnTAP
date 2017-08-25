//
// Created by harta on 8/11/17.
//


//*********************** Includes *****************************
#include <csv.h>
#include <iomanip>
#include "ModRSEM.h"
#include "../ExceptionHandler.h"

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
    _exp_out = (boostFS::path(_expression_outpath)/  _filename).string();
    _rsem_out = _exp_out + ".genes.results";
    if (file_exists(_rsem_out)) {
        print_debug("File found at " + _rsem_out +  "\nmoving to filter transcriptome");
        return std::make_pair(true, "");
    }
    print_debug("File not found at " + _rsem_out +  " Continuing RSEM run.");

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
void ModRSEM::execute(std::map<std::string, QuerySequence> &) {
    // return path
    print_debug("Running RSEM...");

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
    print_debug("Alignment file valid. Preparing reference...");

    ref_exe  = (boostFS::path(_exe_path) / RSEM_PREP_REF_EXE).string();
    ref_path = (boostFS::path(_expression_outpath) / _filename).string() + "_ref";

    // Prepare reference command
    rsem_arg = ref_exe + " "
               + _inpath + " "
               + ref_path;
    std_out = (boostFS::path(_expression_outpath) / _filename).string() + "_rsem_reference";
    print_debug("Executing following command\n" + rsem_arg);
    execute_cmd(rsem_arg.c_str(), std_out);
    print_debug("Reference successfully created");

    print_debug("Running expression analysis...");
    expression_exe = (boostFS::path(_exe_path) / RSEM_CALC_EXP_EXE).string();
    rsem_arg = expression_exe +
               " --" + bam +
               " -p " + std::to_string(_threads) + " " +
               _alignpath +" "+
               ref_path + " " +
               _exp_out;
    if (_ispaired) rsem_arg += " --paired-end";
    std_out = (boostFS::path(_expression_outpath) / _filename).string() + "_rsem_exp";
    print_debug("Executing following command\n" + rsem_arg);
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
std::string ModRSEM::filter(std::map<std::string, QuerySequence> & MAP) {
    print_debug("Beginning to filter transcriptome...");

    unsigned int        count_removed=0;
    unsigned int        count_kept=0;
    unsigned int        count_total=0;      // Used to warn user if high percentage is removed
    float               length;
    float               e_leng;
    float               e_count;
    float               tpm;
    float               fpkm_val;
    float               rejected_percent;
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
    std::stringstream   out_msg;
    GraphingStruct      graphingStruct;

    if (!file_exists(_rsem_out)) {
        throw ExceptionHandler("File does not exist at: " + _rsem_out, ENTAP_ERR::E_RUN_RSEM_EXPRESSION);
    }

    boostFS::remove_all(_processed_path);
    boostFS::create_directories(_processed_path);
    boostFS::create_directories(_figure_path);
    boostFS::path path (_rsem_out);

    fig_txt_box_path = (boostFS::path(_figure_path) / GRAPH_TXT_BOX_PLOT).string();
    fig_png_box_path = (boostFS::path(_figure_path) / GRAPH_PNG_BOX_PLOT).string();

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
    while (in.read_row(geneid, transid, length, e_leng, e_count, tpm, fpkm_val)) {
        count_total++;
        if (fpkm_val > _fpkm) {
            // Kept sequence
            out_file << MAP[geneid].get_sequence() << std::endl;
            QuerySequence *querySequence = &MAP[geneid];
            querySequence->set_is_expression_kept(true);
            querySequence->set_fpkm(fpkm_val);
            file_fig_box << GRAPH_KEPT_FLAG << '\t' <<
                         std::to_string(querySequence->getSeq_length()) << std::endl;
            count_kept++;
        } else {
            // Removed sequence
            removed_file << MAP[geneid].get_sequence() << std::endl;
            file_fig_box << GRAPH_REJECTED_FLAG << '\t'
                         << std::to_string(MAP[geneid].getSeq_length())<<std::endl;
            count_removed++;
        }
    }

    //-----------------------STATISTICS-----------------------//
    out_msg<<std::fixed<<std::setprecision(2);
    out_msg << ENTAP_STATS::SOFTWARE_BREAK
            << "Expression Filtering - RSEM\n"
            << ENTAP_STATS::SOFTWARE_BREAK;
    rejected_percent = ((float)count_removed / count_kept) * 100;
    if (rejected_percent > REJECTED_ERROR_CUTOFF) {
        // Warn user high percentage of transcriptome was rejected
        out_msg << "Warning: A high percentage of the transcriptome was removed: " << rejected_percent;
    }
    out_msg<<"MORE STATS COMING SOON!"<<std::endl;

    out_msg << "Removed: " << count_removed <<
            "\nKept: " << count_kept << std::endl;
    if (count_kept == 0) {
        throw ExceptionHandler("Error in filtering transcriptome, no sequences kept",
                               ENTAP_ERR::E_RUN_RSEM_EXPRESSION);
    }
    out_str = out_msg.str();
    print_statistics(out_str);
    out_file.close();
    removed_file.close();
    file_fig_box.close();
    print_debug("File successfully filtered. Outputs at: " + out_kept + " and: " + out_removed);

    //--------------------------------------------------------//


    //------------------------Graphing------------------------//
    graphingStruct.text_file_path   = fig_txt_box_path;
    graphingStruct.graph_title      = GRAPH_TITLE_BOX_PLOT;
    graphingStruct.fig_out_path     = fig_png_box_path;
    graphingStruct.software_flag    = GRAPH_EXPRESSION_FLAG;
    graphingStruct.graph_type       = GRAPH_BOX_FLAG;
    pGraphingManager->graph(graphingStruct);
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
    print_debug("File is detected to be sam file, running validation");

    std::string rsem_arg;
    std::string out_path;

    rsem_arg = (boostFS::path(_exe_path) / boostFS::path("rsem-sam-validator")).string() +
               " " + _alignpath;
    out_path = (boostFS::path(_expression_outpath) / boostFS::path(filename)).string() + "_rsem_valdate";
    // only thrown in failure in calling rsem
    print_debug("Executing RSEM command:\n" + rsem_arg);
    if (execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) return false;
    print_debug("RSEM validate executed successfully");
    // RSEM does not always return error code if file is invalid, only seen in .err
    return (is_file_empty(out_path+".err"));
}


/**
 * ======================================================================
 * Function ModRSEM::rsem_conv_to_bam(std::string file_name)
 *
 * Description          - Executes convert-sam-for-rsem
 *                      - Converts user inputed file to BAM format
 *
 * Notes                - Not used
 *
 * @param filename      - Transcriptome name just to label output
 *
 * @return              - True if conversion was successful
 *
 * =====================================================================
 */
bool ModRSEM::rsem_conv_to_bam(std::string file_name) {
    print_debug("Converting SAM to BAM");

    std::string bam_out;
    std::string rsem_arg;
    std::string std_out;

    bam_out = (boostFS::path(_expression_outpath) / boostFS::path(file_name)).string();
    rsem_arg = (boostFS::path(_exe_path) / boostFS::path("convert-sam-for-rsem")).string() +
               " -p " + std::to_string(_threads) + " " + _alignpath + " " + bam_out;
    std_out = _expression_outpath + file_name + "_rsem_convert";
    if (execute_cmd(rsem_arg.c_str(), std_out)!=0)return false;
    if (!is_file_empty(std_out+".err")) return false;
    _alignpath = bam_out + ".bam";
    return true;
}


void ModRSEM::set_data(int thread, float fpmk, bool paired) {
    _threads = thread;
    _fpkm = fpmk;
    _ispaired = paired;
}


/**
 * ======================================================================
 * Function ModRSEM::is_file_empty(std::string path)
 *
 * Description          - Check if specific file is empty (has no lines)
 *
 * Notes                - Used for certain RSEM execution that does not
 *                        relay error code of non-zero on failure
 *
 * @param path          - Path to file
 *
 * @return              - True if file is empty
 *
 * =====================================================================
 */
bool ModRSEM::is_file_empty(std::string path) {
    std::ifstream ifstream(path);
    return ifstream.peek() == std::ifstream::traits_type::eof();
}