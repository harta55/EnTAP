//
// Created by harta on 5/7/17.
//

#include <boost/filesystem.hpp>
#include <csv.h>
#include <boost/regex.hpp>
#include <iomanip>
#include "ExpressionAnalysis.h"
#include "ExceptionHandler.h"
#include "EntapInit.h"
#include "EntapConsts.h"
#include "EntapExecute.h"

namespace boostFS = boost::filesystem;

ExpressionAnalysis::ExpressionAnalysis(std::string &input,int t, std::string &exe, std::string &out
    , boost::program_options::variables_map& user_flags) {
    _inpath = input;
    _threads = t;
    _exepath = exe;
    _outpath = out;
    _overwrite = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _ispaired = (bool) user_flags.count("paired-end");
    if (user_flags.count(ENTAP_CONFIG::INPUT_FLAG_ALIGN)) {
        _alignpath = user_flags[ENTAP_CONFIG::INPUT_FLAG_ALIGN].as<std::string>();
    }
    _software_flag = 0;
    _fpkm = user_flags["fpkm"].as<float>();
}

std::string ExpressionAnalysis::execute(std::string input,
                                        std::map<std::string, QuerySequence>& MAP) {
    _inpath = input;
    try {
        switch (_software_flag) {
            case 0:
                return rsem(MAP);
            default:
                return rsem(MAP);
        }
    } catch (const ExceptionHandler &e) {throw e;}
}

std::string ExpressionAnalysis::rsem(std::map<std::string, QuerySequence>& MAP) {
    // return path
    entapInit::print_msg("Running RSEM...");
    boostFS::path out_dir(_outpath+ENTAP_EXECUTE::RSEM_OUT_DIR);
    boostFS::path file_name(_inpath);
    file_name = file_name.stem();
    if (file_name.has_stem()) file_name = file_name.stem(); // for .fasta.fnn
    boostFS::path bam_ext(_alignpath);bam_ext = bam_ext.extension();
    std::string bam = bam_ext.string();
    std::transform(bam.begin(), bam.end(), bam.begin(), ::tolower);
    std::string exp_out_path = out_dir.string() + file_name.string();

    if (_overwrite) {
        boostFS::remove_all(out_dir.c_str());
    } else{
        std::string p = exp_out_path + ".genes.results";
        if (entapInit::file_exists(p)) {
            entapInit::print_msg("File found at " +p +  "\nmoving to filter transcriptome");
            try {
                return rsem_filter(p,MAP);
            } catch (const ExceptionHandler &e) {throw e;}
        }
        entapInit::print_msg("File not found at " + exp_out_path + ".genes.results. Continuing RSEM run.");
    }
    boostFS::create_directories(out_dir);
    if (!entapInit::file_exists(_alignpath)) {
        throw ExceptionHandler("Invalid file path for BAM/SAM file, exiting...",
                               ENTAP_ERR::E_INIT_TAX_READ);
    }
    std::string rsem_arg,out_path;
    if (bam.compare(".sam")==0) {
        entapInit::print_msg("File is detected to be sam file, running validation "
                                     "and conversion to bam");
        rsem_arg = _exepath + "rsem-sam-validator " + _alignpath;
        out_path = out_dir.string() + file_name.string() + "_rsem_valdate";
        if (entapInit::execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) {
            // only thrown in failure in calling rsem
            throw ExceptionHandler("Error in validating sam file", ENTAP_ERR::E_INIT_TAX_READ);
        }
        // RSEM does not return error code if file is invalid, only seen in .err
        if (!is_file_empty(out_path+".err")) {
            throw ExceptionHandler("Alignment file invalid!", ENTAP_ERR::E_INIT_TAX_READ);
        }
        entapInit::print_msg("Alignment file valid. Converting to BAM");
        std::string bam_out = out_dir.string() + file_name.string();
        rsem_arg = _exepath+ "convert-sam-for-rsem " + " -p " + std::to_string(_threads)
                   + " "+_alignpath + " " + bam_out;
        out_path = out_dir.string() + file_name.string() + "_rsem_convert";
        if (entapInit::execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) {
            // execution error, dif from conversion error
            throw ExceptionHandler("Error in converting sam file", ENTAP_ERR::E_INIT_TAX_READ);
        }
        if (!is_file_empty(out_path+".err")) {
            throw ExceptionHandler("Error in converting sam file", ENTAP_ERR::E_INIT_TAX_READ);
        }
        _alignpath = bam_out + ".bam";

    } else if(bam.compare(".bam")==0) {
        entapInit::print_msg("File is detected to be bam file, validating...");
        rsem_arg = _exepath + "rsem-sam-validator " + _alignpath;
        out_path = out_dir.string() + file_name.string() + "_rsem_valdate";
        if (entapInit::execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) {
            throw ExceptionHandler("Error in validating bam file", ENTAP_ERR::E_INIT_TAX_READ);
        }
        if (!is_file_empty(out_path+".err")) {
            throw ExceptionHandler("Alignment file invalid!", ENTAP_ERR::E_INIT_TAX_READ);
        }
        entapInit::print_msg("Alignment file valid. Continuing...");
    } else {
        throw ExceptionHandler("Unknown extension found in the alignment file",
                               ENTAP_ERR::E_INIT_TAX_READ);
    }
    // Now have valid BAM file to run rsem
    entapInit::print_msg("Preparing reference");
    std::string ref_out_path = out_dir.string() + file_name.string() + "_ref";
    boostFS::create_directories(ref_out_path);
    rsem_arg = _exepath + "rsem-prepare-reference " + _inpath +
               " " + ref_out_path+"/"+file_name.string();
    std::string std_out = out_dir.string() + file_name.string() + "_rsem_reference";
    entapInit::print_msg("Executing following command\n" + rsem_arg);
    entapInit::execute_cmd(rsem_arg.c_str(), std_out);
    entapInit::print_msg("Reference successfully created");

    entapInit::print_msg("Running expression analysis...");
    rsem_arg = _exepath + "rsem-calculate-expression " +
               "--bam " + "-p " + std::to_string(_threads) + " " + _alignpath +" "+ ref_out_path +
               "/"+file_name.string()+ " " +exp_out_path;
    if (_ispaired) rsem_arg += " --paired-end";
    std_out = out_dir.string() + file_name.string() + "_rsem_exp";
    entapInit::print_msg("Executing following command\n" + rsem_arg);
    if (entapInit::execute_cmd(rsem_arg.c_str(), std_out)!=0) {
        throw ExceptionHandler("Error in running expression analysis",ENTAP_ERR::E_INIT_TAX_READ);
    }
    std::string out_results = exp_out_path + ".genes.results";
    try {
        return rsem_filter(out_results,MAP);
    } catch (const ExceptionHandler &e) {throw e;}
}

std::string ExpressionAnalysis::rsem_filter(std::string &results_path,
                                            std::map<std::string, QuerySequence>& MAP) {
    entapInit::print_msg("Beginning to filter transcriptome...");
    if (!entapInit::file_exists(results_path)) {
        throw ExceptionHandler("File does not exist at: " + results_path,
            ENTAP_ERR::E_RUN_RSEM_EXPRESSION);
    }
    std::stringstream out_msg;out_msg<<std::fixed<<std::setprecision(2);
    out_msg << ENTAP_STATS::SOFTWARE_BREAK
            << "Expression Filtering - RSEM\n"
            << ENTAP_STATS::SOFTWARE_BREAK;
    out_msg<<"MORE STATES COMING SOON!"<<std::endl;

    std::string process_dir = _outpath + ENTAP_EXECUTE::RSEM_OUT_DIR + "processed";
    boostFS::remove_all(process_dir); boostFS::create_directories(process_dir);
    boostFS::path path (results_path);
    std::string out_path = _outpath +
                           path.stem().stem().string() + "_removed.fasta";
    std::string out_removed = _outpath +
                              path.stem().string() + "_filtered.fasta";
    std::ofstream out_file(out_path, std::ios::out | std::ios::app);
    std::ofstream removed_file(out_removed, std::ios::out | std::ios::app);
    io::CSVReader<ENTAP_EXECUTE::RSEM_COL_NUM, io::trim_chars<' '>,
    io::no_quote_escape<'\t'>> in(results_path);
    in.next_line();
    long count_removed=0, count_kept=0;
    std::string geneid, transid;
    float length, e_leng, e_count, tpm, fpkm_val;
    while (in.read_row(geneid, transid, length, e_leng, e_count, tpm, fpkm_val)) {
        if (fpkm_val > _fpkm) {
            out_file << MAP[geneid].get_sequence() << std::endl;
            MAP[geneid].set_is_expression_kept(true);
            count_kept++;
        } else {
            removed_file << MAP[geneid].get_sequence() << std::endl;
            count_removed++;
        }
    }

    out_msg << "Removed: " << count_removed <<
            "\nKept: " << count_kept << std::endl;
    if (count_kept == 0) {
        throw ExceptionHandler("Error in filtering transcriptome, no sequences kept",
        ENTAP_ERR::E_RUN_RSEM_EXPRESSION);
    }
    std::string out_str = out_msg.str();
    entapExecute::print_statistics(out_str,_outpath);
    out_file.close();removed_file.close();
    entapInit::print_msg("File successfully filtered. Outputs at: " + out_path + " and: " +
                         out_removed);
    return out_path;
}

bool ExpressionAnalysis::is_file_empty(std::string path) {
    std::ifstream ifstream(path);
    return ifstream.peek() == std::ifstream::traits_type::eof();
}

