//
// Created by harta on 5/10/17.
//

#include <boost/filesystem.hpp>
#include "SimilaritySearch.h"
#include "ExceptionHandler.h"
#include "EntapInit.h"
#include "EntapConsts.h"

namespace boostFS = boost::filesystem;

SimilaritySearch::SimilaritySearch(std::list<std::string> &databases, std::string input,
                           int threads, bool overwrite, std::string exe, std::string out) {
    this->_database_paths = databases;
    this->_input = input;
    this->_threads = threads;
    this->_overwrite = overwrite;
    this->_exe = exe;
    this->_outpath = out;
}

std::list<std::string> SimilaritySearch::execute(short software, std::string updated_input,
    bool blast) {
    this->_input = updated_input;
    this->_blastp = blast;
    try {
        switch (software) {
            case 0:
                return diamond();
            default:
                return diamond();
        }
    } catch (ExceptionHandler &e) {throw e;}
    return std::list<std::string>();
}

std::list<std::string> SimilaritySearch::diamond() {
    // not always known (depending on starting state)
    entapInit::print_msg("Beginning to execute DIAMOND...");
    std::string diamond_out = _outpath + ENTAP_CONFIG::DIAMOND_DIR;
    std::string blast_type;
    _blastp ? blast_type = "blastp" : blast_type = "blastx";
    std::list<std::string> out_paths;
    if (!entapInit::file_exists(_input)) {
        throw ExceptionHandler("Transcriptome file not found",ENTAP_ERR::E_INIT_INDX_DATA_NOT_FOUND);
    }
    boostFS::path transc_name(_input); transc_name=transc_name.stem();
    if (transc_name.has_stem()) transc_name = transc_name.stem(); //.fasta.faa

    if (_overwrite) {
        boostFS::remove_all(diamond_out.c_str());
    } else {
        std::list<std::string> temp = verify_diamond_files(diamond_out,blast_type,
                                                           transc_name.string());
        if (!temp.empty()) return temp;
    }
    boostFS::create_directories(diamond_out);

    // database verification already ran, don't need to verify each path
    try {
        for (std::string data_path : _database_paths) {
            // assume all paths should be .dmnd
            boostFS::path file_name(data_path); file_name=file_name.stem();
            entapInit::print_msg("Searching against database located at: " +
                                 data_path + "...");
            std::string out_path = diamond_out + blast_type + "_" + transc_name.string() + "_" +
                                   file_name.string() + ".out";
            std::string std_out = diamond_out + blast_type + "_" + transc_name.string() + "_" +
                                  file_name.string() + "_std";
            diamond_blast(_input, out_path, std_out,data_path, _threads, blast_type);
            entapInit::print_msg("Success! Results written to " + out_path);
            out_paths.push_back(out_path);
        }
    } catch (ExceptionHandler &e) {
        throw ExceptionHandler(e.what(), e.getErr_code());
    }
    return out_paths;
}

void SimilaritySearch::diamond_blast(std::string input_file, std::string output_file, std::string std_out,
                   std::string &database,int &threads, std::string &blast) {
    std::string diamond_run = _exe + " " + blast +" -d " + database +
                              " --sensitive" + " -k 10" + " -q " + input_file + " -o " + output_file + " -p " + std::to_string(threads) +" -f " +
                              "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle";
    entapInit::print_msg("\nExecuting Diamond:\n" + diamond_run);
    if (entapInit::execute_cmd(diamond_run, std_out) != 0) {
        throw ExceptionHandler("Error in DIAMOND run with database located at: " +
                               database, ENTAP_ERR::E_INIT_TAX_INDEX);
    }
}

std::list<std::string> SimilaritySearch::verify_diamond_files(std::string &outpath,
                std::string &blast, std::string name) {
    entapInit::print_msg("Override unselected, checking for diamond files"
                                 " of selected databases...");
    std::list<std::string> out_list;
    for (std::string data_path : _database_paths) {
        // assume all paths should be .dmnd
        boostFS::path file_name(data_path);
        file_name = file_name.stem();
        std::string temp_out = outpath + blast + "_" + name + "_" +
                               file_name.string() + ".out";
        if (!entapInit::file_exists(temp_out)){
            entapInit::print_msg("File at: " + temp_out + " not found, running diamond");
            out_list.clear();return out_list;
        }
        out_list.push_back(temp_out);
    }
    entapInit::print_msg("All diamond files found, skipping this stage of enTAP");
    return out_list;
}
