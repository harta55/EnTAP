//
// Created by harta on 5/7/17.
//

#include <boost/filesystem.hpp>
#include <csv.h>
#include <boost/regex.hpp>
#include <iomanip>
#include "ExpressionAnalysis.h"
#include "ExceptionHandler.h"
#include "EntapConfig.h"
#include "EntapGlobals.h"
#include "EntapExecute.h"

namespace boostFS = boost::filesystem;

ExpressionAnalysis::ExpressionAnalysis(std::string &input,int t, std::string &exe, std::string &out
    , boost::program_options::variables_map& user_flags) {
    _inpath = input;
    _threads = t;
    _exepath = exe;
    _outpath = out;
    _software_flag = 0;
    _overwrite = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _ispaired = (bool) user_flags.count("paired-end");
    _alignpath = user_flags[ENTAP_CONFIG::INPUT_FLAG_ALIGN].as<std::string>(); // Checked in EntapExecute.cpp
    _fpkm = user_flags[ENTAP_CONFIG::INPUT_FLAG_FPKM].as<float>();
    _rsem_dir = (boostFS::path(out) / boostFS::path(RSEM_OUT_DIR)).string();
    _proc_dir = (boostFS::path(_rsem_dir) / boostFS::path(RSEM_PROCESSED_DIR)).string();
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
    print_debug("Running RSEM...");

    boostFS::path                   file_name(_inpath);
    boostFS::path                   bam_ext(_alignpath);
    std::string                     bam;
    std::string                     rsem_arg;
    std::string                     out_path;

    while (file_name.has_extension()) file_name = file_name.stem();
    bam_ext = bam_ext.extension();
    bam = bam_ext.string();
    std::transform(bam.begin(), bam.end(), bam.begin(), ::tolower);
    std::string exp_out_path = _rsem_dir + file_name.string();

    if (_overwrite) {
        boostFS::remove_all(_rsem_dir);
    } else{
        std::string p = exp_out_path + ".genes.results";
        if (file_exists(p)) {
            print_debug("File found at " +p +  "\nmoving to filter transcriptome");
            try {
                return rsem_filter(p,MAP);
            } catch (const ExceptionHandler &e) {throw e;}
        }
        print_debug("File not found at " + exp_out_path + ".genes.results. Continuing RSEM run.");
    }
    boostFS::create_directories(_rsem_dir);

    if (!file_exists(_alignpath)) {
        throw ExceptionHandler("Invalid file path for BAM/SAM file, exiting...", ENTAP_ERR::E_INIT_TAX_READ);
    }
    //todo separate into methods
    if (bam.compare(".sam")==0) {
        if (!rsem_validate_file(file_name.string())){
            throw ExceptionHandler("SAM file invalid", ENTAP_ERR::E_RUN_RSEM_VALIDATE);
        }
        if (!rsem_conv_to_bam(file_name.string())) {
            throw ExceptionHandler("Error in converting sam file", ENTAP_ERR::E_RUN_RSEM_CONVERT);
        }

    } else if(bam.compare(".bam")==0) {
        print_debug("File is detected to be bam file, validating...");
        rsem_arg = _exepath + "rsem-sam-validator " + _alignpath;
        out_path = _rsem_dir + file_name.string() + "_rsem_valdate";
        if (execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) {
            throw ExceptionHandler("Error in validating bam file", ENTAP_ERR::E_INIT_TAX_READ);
        }
        if (!is_file_empty(out_path+".err")) {
            throw ExceptionHandler("Alignment file invalid!", ENTAP_ERR::E_INIT_TAX_READ);
        }
        print_debug("Alignment file valid. Continuing...");
    } else {
        throw ExceptionHandler("Unknown extension found in the alignment file",
                               ENTAP_ERR::E_INIT_TAX_READ);
    }

    // Now have valid BAM file to run rsem
    print_debug("Preparing reference");
    std::string ref_out_path = _rsem_dir + file_name.string() + "_ref";
    boostFS::create_directories(ref_out_path);
    rsem_arg = _exepath + "rsem-prepare-reference " + _inpath +
               " " + ref_out_path+"/"+file_name.string();
    std::string std_out = _rsem_dir + file_name.string() + "_rsem_reference";
    print_debug("Executing following command\n" + rsem_arg);
    execute_cmd(rsem_arg.c_str(), std_out);
    print_debug("Reference successfully created");

    print_debug("Running expression analysis...");
    rsem_arg = _exepath + "rsem-calculate-expression " +
               "--bam " + "-p " + std::to_string(_threads) + " " + _alignpath +" "+ ref_out_path +
               "/"+file_name.string()+ " " +exp_out_path;
    if (_ispaired) rsem_arg += " --paired-end";
    std_out = _rsem_dir + file_name.string() + "_rsem_exp";
    print_debug("Executing following command\n" + rsem_arg);
    if (execute_cmd(rsem_arg.c_str(), std_out)!=0) {
        throw ExceptionHandler("Error in running expression analysis",ENTAP_ERR::E_INIT_TAX_READ);
    }
    std::string out_results = exp_out_path + ".genes.results";
    try {
        return rsem_filter(out_results,MAP);
    } catch (const ExceptionHandler &e) {throw e;}
}

std::string ExpressionAnalysis::rsem_filter(std::string &results_path,
                                            std::map<std::string, QuerySequence>& MAP) {
    print_debug("Beginning to filter transcriptome...");

    long                count_removed=0;
    long                count_kept=0;
    float               length;
    float               e_leng;
    float               e_count;
    float               tpm;
    float               fpkm_val;
    std::string         geneid;
    std::string         transid;
    std::string         out_str;
    std::stringstream   out_msg;

    if (!file_exists(results_path)) {
        throw ExceptionHandler("File does not exist at: " + results_path,
            ENTAP_ERR::E_RUN_RSEM_EXPRESSION);
    }

    boostFS::remove_all(_proc_dir);
    boostFS::create_directories(_proc_dir);
    boostFS::path path (results_path);
    std::string out_path = _outpath +
                           path.stem().stem().string() + "_removed.fasta";
    std::string out_removed = _outpath +
                              path.stem().string() + "_filtered.fasta";
    std::ofstream out_file(out_path, std::ios::out | std::ios::app);
    std::ofstream removed_file(out_removed, std::ios::out | std::ios::app);
    io::CSVReader<RSEM_COL_NUM, io::trim_chars<' '>,
    io::no_quote_escape<'\t'>> in(results_path);
    in.next_line();
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

    //------------------STATISTICS--------------------//
    out_msg<<std::fixed<<std::setprecision(2);
    out_msg << ENTAP_STATS::SOFTWARE_BREAK
            << "Expression Filtering - RSEM\n"
            << ENTAP_STATS::SOFTWARE_BREAK;
    out_msg<<"MORE STATS COMING SOON!"<<std::endl;

    out_msg << "Removed: " << count_removed <<
            "\nKept: " << count_kept << std::endl;
    if (count_kept == 0) {
        throw ExceptionHandler("Error in filtering transcriptome, no sequences kept",
        ENTAP_ERR::E_RUN_RSEM_EXPRESSION);
    }
    out_str = out_msg.str();
    print_statistics(out_str);
    out_file.close();removed_file.close();
    print_debug("File successfully filtered. Outputs at: " + out_path + " and: " +
                         out_removed);
    return out_path;
}

bool ExpressionAnalysis::is_file_empty(std::string path) {
    std::ifstream ifstream(path);
    return ifstream.peek() == std::ifstream::traits_type::eof();
}

bool ExpressionAnalysis::rsem_validate_file(std::string filename) {
    print_debug("File is detected to be sam file, running validation");

    std::string rsem_arg;
    std::string out_path;

    rsem_arg = (boostFS::path(_exepath) / boostFS::path("rsem-sam-validator")).string() +
            " " + _alignpath;
    out_path = (boostFS::path(_rsem_dir) / boostFS::path(filename)).string() + "_rsem_valdate";
    // only thrown in failure in calling rsem
    if (execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) return false;
    // RSEM does not return error code if file is invalid, only seen in .err
    return (!is_file_empty(out_path+".err"));
}

bool ExpressionAnalysis::rsem_conv_to_bam(std::string file_name) {
    print_debug("Converting SAM to BAM");

    std::string bam_out;
    std::string rsem_arg;
    std::string std_out;

    bam_out = (boostFS::path(_rsem_dir) / boostFS::path(file_name)).string();
    rsem_arg = (boostFS::path(_exepath) / boostFS::path("convert-sam-for-rsem")).string() +
                " -p " + std::to_string(_threads) + " " + _alignpath + " " + bam_out;
    std_out = _rsem_dir + file_name + "_rsem_convert";
    if (execute_cmd(rsem_arg.c_str(), std_out)!=0)return false;
    if (!is_file_empty(std_out+".err")) return false;
    _alignpath = bam_out + ".bam";
    return true;
}

