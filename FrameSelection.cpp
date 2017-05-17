//
// Created by harta on 5/7/17.
//

#include <boost/filesystem/operations.hpp>
#include <fstream>
#include "FrameSelection.h"
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapInit.h"
#include "EntapConsts.h"

namespace boostFS = boost::filesystem;

// can accept version/other for dependency injection
FrameSelection::FrameSelection(std::string &input, std::string &exe, std::string &out, bool overwrite) {
    this->_exe_path = exe;
    this->_outpath = out;
    this->_inpath = input;
    this->_overwrite = overwrite;
}

std::string FrameSelection::execute(short software) {
    try {
        switch (software) {
            case 0:
                return genemarkst();
            default:
                return genemarkst();
        }
    } catch (ExceptionHandler &e) {throw e;}
}

std::string FrameSelection::genemarkst() {
    // Outfiles: file/path.faa, file/path.fnn
    // assumes working directory as output right now
    boost::filesystem::path file_name(_inpath); file_name = file_name.filename();
    std::list<std::string> out_names {file_name.string()+".faa",
                                      file_name.string()+".fnn"};
    std::string genemark_out_dir = _outpath + "genemark/";
    std::string final_out = genemark_out_dir + file_name.string() + ".faa";

    if (_overwrite) {
        boostFS::remove_all(genemark_out_dir.c_str());
    } else {
        if (entapInit::file_exists(final_out)) {
            entapInit::print_msg("File found at: " + final_out + "\n"
                 "continuing enTAP with this file and skipping frame selection");
            return final_out;
        }
        entapInit::print_msg("File not found at " + final_out + " so continuing frame selection");
    }
    boostFS::create_directories(genemark_out_dir);

    std::string genemark_cmd = _exe_path + " -faa -fnn " + _inpath;
    std::string genemark_std_out = genemark_out_dir + "genemark_run";
    entapInit::print_msg("Running genemark...\n" + genemark_cmd);

    if (entapInit::execute_cmd(genemark_cmd,genemark_std_out) != 0 ) {
        throw ExceptionHandler("Error in running genemark at file located at: " +
                               _inpath, ENTAP_ERR::E_INIT_INDX_DATA_NOT_FOUND);
    }
    entapInit::print_msg("Success!");

    // Format genemarks-t output (remove blank lines)
    entapInit::print_msg("Formatting genemark files");

    std::string line;
    for (std::string path : out_names) {
        std::ifstream in_file(path);
        std::string temp_name = path+"_alt";
        std::string out_path = genemark_out_dir+path;
        std::ofstream out_file(path+"_alt");
        while (getline(in_file,line)){
            if (!line.empty()) {
                out_file << line << '\n';
            }
        }
        in_file.close();
        out_file.close();
        if (remove(path.c_str())!=0 || rename(temp_name.c_str(),out_path.c_str())!=0) {
            throw ExceptionHandler("Error formatting/moving genemark results", ENTAP_ERR::E_INIT_TAX_READ);
        }
    }
    std::string lst_file = file_name.string() + ".lst";
    std::string out_lst = genemark_out_dir + lst_file;
    std::string out_gmst_log = genemark_out_dir+ENTAP_EXECUTE::GENEMARK_LOG_FILE;

    if (rename(lst_file.c_str(),out_lst.c_str())!=0 ||
            rename(ENTAP_EXECUTE::GENEMARK_LOG_FILE.c_str(),out_gmst_log.c_str())!=0) {
        throw ExceptionHandler("Error moving genemark results", ENTAP_ERR::E_INIT_TAX_READ);
    }
    entapInit::print_msg("Success!");
    genemarkStats(final_out);
    return final_out;
}

/**
 *
 * @param protein_path - path to .faa file
 * @return
 */
std::string FrameSelection::genemarkStats(std::string &protein_path) {
    // generate maps, query->sequence
    std::map<std::string,std::string> protein_map = entapExecute::generate_seq_map(protein_path);
    std::string out_removed_path = _outpath + ENTAP_EXECUTE::GENEMARK_OUT_PATH +
        ENTAP_EXECUTE::GENEMARK_OUT_UNSELECTED;
    float min_removed,min_selected,max_removed,max_selected,avg_removed,avg_selected;
    unsigned int count_total = 0, count_removed = 0;

    struct frame_seq {
        unsigned int length;
        std::string sequence;
        std::string frame_type;
    };
    std::map<std::string,frame_seq> nucleotide_map;

    std::ifstream in_file(_inpath);
    std::ofstream removed_file(out_removed_path, std::ios::out | std::ios::app);
    std::string line, sequence, seq_id;
    frame_seq nucleotide_sequence;
    while (getline(in_file,line)){
        if (line.empty()) continue;
        if (line.find(">")==0) {
            count_total++;
            if (!seq_id.empty()) {
                nucleotide_map.emplace(seq_id,nucleotide_sequence);
                if (protein_map.find(seq_id) == protein_map.end()) {
                    // Frame was not found
                    count_removed++;
                    removed_file<<nucleotide_sequence.sequence;
                }
            }
            seq_id = line.substr(line.find(">")+1,line.find(" "));
            nucleotide_sequence.sequence = line + "\n";
        } else {
            nucleotide_sequence.sequence += line + "\n";
        }
    }
    in_file.close(); removed_file.close();
    return "";
}

