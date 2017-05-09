//
// Created by harta on 5/7/17.
//

#include <boost/filesystem.hpp>
#include "ExpressionAnalysis.h"
#include "ExceptionHandler.h"
#include "EntapInit.h"
#include "EntapConsts.h"
#include "EntapExecute.h"

namespace boostFS = boost::filesystem;

ExpressionAnalysis::ExpressionAnalysis(std::string &input,int t, std::string &exe, std::string &out) {
    this->_inpath = input;
//    this->_alignpath = alignment;
//    this->_ispaired = paired;
    this->_threads = t;
    this->_exepath = exe;
    this->_outpath = out;
}

std::string ExpressionAnalysis::execute(short software, bool paired, std::string align) {
    this->_alignpath = align;
    this->_ispaired = paired;
    try {
        switch (software) {
            case 0:
                return rsem();
            default:
                return rsem();
        }
    } catch (ExceptionHandler &e) {throw e;}
}

std::string ExpressionAnalysis::rsem() {
    // return path
    entapInit::print_msg("Running RSEM...");
    boostFS::path out_dir(_outpath+"rsem/");
    boostFS::remove_all(out_dir.c_str());
    boostFS::create_directories(out_dir);
    boostFS::path file_name(_inpath);
    file_name = file_name.stem();
    if (file_name.has_stem()) file_name = file_name.stem(); // for .fasta.fnn
    boostFS::path bam_ext(_alignpath);bam_ext = bam_ext.extension();
    std::string bam = bam_ext.string();
    std::transform(bam.begin(), bam.end(), bam.begin(), ::tolower);

    if (_alignpath.empty()) {
        entapInit::print_msg("No BAM/SAM file provided, exiting RSEM run");
        return _inpath;
    }
    if (!entapInit::file_exists(_alignpath)) {
        throw ExceptionHandler("Invalid file path for BAM/SAM file, exiting...",
                               ENTAP_ERR::E_INIT_TAX_READ);
    }
    std::string rsem_arg;
    std::string out_path;
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
        if (!entapExecute::is_file_empty(out_path+".err")) {
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
        if (!entapExecute::is_file_empty(out_path+".err")) {
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
        if (!entapExecute::is_file_empty(out_path+".err")) {
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
    std::string exp_out_path = out_dir.string() + file_name.string();
    rsem_arg = _exepath + "rsem-calculate-expression " + "--paired-end " +
               "--bam " + "-p " + std::to_string(_threads) + " " + _alignpath +" "+ ref_out_path +
               "/"+file_name.string()+ " " +exp_out_path;
    if (_ispaired) rsem_arg += " --paired-end";
    std_out = out_dir.string() + file_name.string() + "_rsem_exp";
    entapInit::print_msg("Executing following command\n" + rsem_arg);
    if (entapInit::execute_cmd(rsem_arg.c_str(), std_out)!=0) {
        throw ExceptionHandler("Error in running expression analysis",ENTAP_ERR::E_INIT_TAX_READ);
    }
    return out_path + ".genes.results";
}

