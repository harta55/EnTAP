//
// Created by harta on 5/7/17.
//

#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <boost/regex.hpp>
#include <iomanip>
#include "FrameSelection.h"
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapInit.h"
#include "EntapConsts.h"

namespace boostFS = boost::filesystem;

// can accept version/other for dependency injection
FrameSelection::FrameSelection(std::string &input, std::string &exe, std::string &out,
                               boost::program_options::variables_map &user_flags) {
    _exe_path = exe;
    _outpath = out;
    _inpath = input;
    _overwrite = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _software_flag = 0;
}

std::string FrameSelection::execute(std::map<std::string,QuerySequence> &SEQUENCES) {
    try {
        switch (_software_flag) {
            case 0:
                return genemarkst(SEQUENCES);
            default:
                return genemarkst(SEQUENCES);
        }
    } catch (ExceptionHandler &e) {throw e;}
}

std::string FrameSelection::genemarkst(std::map<std::string,QuerySequence> &SEQUENCES) {
    // Outfiles: file/path.faa, file/path.fnn
    // assumes working directory as output right now
    boost::filesystem::path file_name(_inpath); file_name = file_name.filename();
    std::list<std::string> out_names {file_name.string()+".faa",
                                      file_name.string()+".fnn"};
    std::string genemark_out_dir = _outpath + ENTAP_EXECUTE::GENEMARK_OUT_PATH;
    std::string final_out = genemark_out_dir + file_name.string() + ".faa";
    std::string lst_file = file_name.string() + ".lst";
    std::string out_lst = genemark_out_dir + lst_file;
    std::string out_gmst_log = genemark_out_dir+ENTAP_EXECUTE::GENEMARK_LOG_FILE;
    std::string out_hmm_file = genemark_out_dir+ENTAP_EXECUTE::GENEMARK_HMM_FILE;

    if (_overwrite) {
        boostFS::remove_all(genemark_out_dir.c_str());
    } else {
        if (entapInit::file_exists(final_out) && entapInit::file_exists(out_lst)) {
            entapInit::print_msg("File found at: " + final_out + "\n"
                 "continuing enTAP with this file and skipping frame selection");
            try {
                genemarkStats(final_out,out_lst,SEQUENCES);
            } catch (ExceptionHandler &e){throw e;}
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
        in_file.close();out_file.close();
        if (remove(path.c_str())!=0 || rename(temp_name.c_str(),out_path.c_str())!=0) {
            throw ExceptionHandler("Error formatting/moving genemark results", ENTAP_ERR::E_INIT_TAX_READ);
        }
    }
    if (rename(lst_file.c_str(),out_lst.c_str())!=0 ||
        rename(ENTAP_EXECUTE::GENEMARK_LOG_FILE.c_str(),out_gmst_log.c_str())!=0 ) {
        throw ExceptionHandler("Error moving genemark results", ENTAP_ERR::E_INIT_TAX_READ);
    }
    if (entapInit::file_exists(ENTAP_EXECUTE::GENEMARK_HMM_FILE)) {
        rename(ENTAP_EXECUTE::GENEMARK_HMM_FILE.c_str(),out_hmm_file.c_str());
    }
    entapInit::print_msg("Success!");
    try {
        genemarkStats(final_out,out_lst,SEQUENCES);
    } catch (ExceptionHandler &e){throw e;}
    return final_out;
}

/**
 *
 * @param protein_path - path to .faa file
 * @param lst_path - path to .lst genemark file
 */
void FrameSelection::genemarkStats(std::string &protein_path, std::string &lst_path,
                                   std::map<std::string,QuerySequence> &SEQUENCES) {
    // generate maps, query->sequence
    entapInit::print_msg("Beginning to calculate Genemark statistics...");
    std::string processed_path = _outpath + ENTAP_EXECUTE::FRAME_SELECTION_PROCESSED +"/";
    boostFS::create_directories(processed_path);
    std::string out_removed_path = processed_path + ENTAP_EXECUTE::FRAME_SELECTION_LOST;
    std::string out_internal_path = processed_path + ENTAP_EXECUTE::FRAME_SELECTION_INTERNAL;
    std::string out_complete_path = processed_path + ENTAP_EXECUTE::FRAME_SELECTION_COMPLTE;
    std::string out_partial_path = processed_path + ENTAP_EXECUTE::FRAME_SELECTION_PARTIAL;

    // all nucleotide lengths
    unsigned long min_removed=10000,min_selected=10000,max_removed=0,max_selected=0,
            total_removed_len=0,total_kept_len=0, count_selected = 0, count_removed = 0,
            count_partial_5=0, count_partial_3=0, count_internal=0, count_complete=0;
    std::string min_removed_seq, min_kept_seq, max_removed_seq, max_kept_seq;
    try {
        std::map<std::string,frame_seq> protein_map = genemark_parse_protein(protein_path);
        genemark_parse_lst(lst_path,protein_map);

        std::string lost_flag = "lost";
        std::string selected_flag = "selected";
        std::map<std::string, std::ofstream*> file_map;
        file_map[lost_flag] =
                new std::ofstream(out_removed_path, std::ios::out | std::ios::app);
        file_map[ENTAP_EXECUTE::FRAME_SELECTION_INTERNAL_FLAG] =
                new std::ofstream(out_internal_path, std::ios::out | std::ios::app);
        file_map[ENTAP_EXECUTE::FRAME_SELECTION_COMPLETE_FLAG] =
                new std::ofstream(out_complete_path, std::ios::out | std::ios::app);
        file_map[ENTAP_EXECUTE::FRAME_SELECTION_FIVE_FLAG] =
                new std::ofstream(out_partial_path, std::ios::out | std::ios::app);
        file_map[ENTAP_EXECUTE::FRAME_SELECTION_THREE_FLAG] = file_map[ENTAP_EXECUTE::FRAME_SELECTION_FIVE_FLAG];

        std::map<std::string, unsigned long> count_map ={
                {ENTAP_EXECUTE::FRAME_SELECTION_INTERNAL_FLAG,count_internal},
                {ENTAP_EXECUTE::FRAME_SELECTION_COMPLETE_FLAG,count_complete},
                {ENTAP_EXECUTE::FRAME_SELECTION_FIVE_FLAG,count_partial_5},
                {ENTAP_EXECUTE::FRAME_SELECTION_THREE_FLAG,count_partial_3},
        };

        std::vector<unsigned long> all_kept_lengths, all_lost_lengths;
        for (auto& pair : SEQUENCES) {
            std::map<std::string,frame_seq>::iterator p_it = protein_map.find(pair.first);
            if (p_it != protein_map.end()) {
                count_selected++;
                pair.second.setSequence(p_it->second.sequence);
                pair.second.setFrame(p_it->second.frame_type);

                std::string sequence = p_it->second.sequence;
                std::string frame_type = p_it->second.frame_type;
                unsigned long length = pair.second.getSeq_length();  // Nucleotide sequence length

                if (length < min_selected) {
                    min_selected = length;
                    min_kept_seq = pair.first;
                }
                if (length > max_selected) {
                    max_selected = length;
                    max_kept_seq = pair.first;
                }
                total_kept_len += length;
                all_kept_lengths.push_back(length);
                std::map<std::string, std::ofstream*>::iterator file_it = file_map.find(frame_type);
                if (file_it != file_map.end()) {
                    *file_it->second << sequence;
                    count_map[frame_type]++;
                } else {
                    throw ExceptionHandler("Unknown frame flag found", ENTAP_ERR::E_RUN_GENEMARK_STATS);
                }

            } else {
                // Lost sequence
                count_removed++;
                *file_map[lost_flag] << pair.second.getSequence() ;
                unsigned long length = pair.second.getSeq_length();  // Nucleotide sequence length

                if (length < min_removed) {
                    min_removed = length;
                    min_removed_seq = pair.first;
                }
                if (length > max_removed) {
                    max_removed_seq = pair.first;
                    max_removed = length;
                }
                all_lost_lengths.push_back(length);
                total_removed_len += length;
            }
        }
        for(auto& pair : file_map) {
            if (pair.first.compare(ENTAP_EXECUTE::FRAME_SELECTION_THREE_FLAG)!=0){
                pair.second->close();
                delete pair.second;
                pair.second = 0;
            }
        }
        // Calculate and print stats
        double avg_selected = (double)total_kept_len / count_selected;
        double avg_lost = (double)total_removed_len / count_removed;
        std::stringstream stat_output; stat_output<<std::fixed<<std::setprecision(2);
        stat_output <<
                    ENTAP_STATS::SOFTWARE_BREAK <<
                    "Frame Selection: GenemarkS-T" <<
                    "\n"<<ENTAP_STATS::SOFTWARE_BREAK;
        stat_output <<
                    "Total sequences frame selected: "                 << count_selected            <<
                    "\n\tThese protein sequences were written to: "    << protein_path              <<
                    "\nTotal sequences lost during selection: "        << count_removed             <<
                    "\n\tThese nucleotide sequences were written to: " << out_removed_path          <<
                    "\nThere were " + count_map[ENTAP_EXECUTE::FRAME_SELECTION_FIVE_FLAG]           <<
                    " 5 prime partials and " + count_map[ENTAP_EXECUTE::FRAME_SELECTION_THREE_FLAG] <<
                    " 3 prime partials" <<
                    "\n\tAll partials were written to: " + out_partial_path <<
                    "\nThere were " + count_map[ENTAP_EXECUTE::FRAME_SELECTION_COMPLETE_FLAG]       <<
                    " complete genes\n\tThese were written to: " << out_complete_path               <<
                    "\nThere were " << count_map[ENTAP_EXECUTE::FRAME_SELECTION_INTERNAL_FLAG]      <<
                    " internal genes\n\tThese were written to: " << out_internal_path               <<"\n\n";

        std::pair<unsigned long, unsigned long> kept_n =
                entapExecute::calculate_N_vals(all_kept_lengths,total_kept_len);
        stat_output << "New transcriptome reference:"
                    "\n\tTotal sequences: "      << count_selected <<
                    "\n\tTotal length of transcriptome(bp): "      << total_kept_len <<
                    "\n\tAverage length(bp): "   << avg_selected   <<
                    "\n\tn50: "                  << kept_n.first   <<
                    "\n\tn90: "                  << kept_n.second  <<
                    "\n\tLongest sequence(bp): " << max_selected   << " (" << max_kept_seq << ")" <<
                    "\n\tShortest sequence(bp): "<< min_selected   << " (" << min_kept_seq << ")";

        if (count_removed > 0) {
            std::pair<unsigned long, unsigned long> removed_n =
                    entapExecute::calculate_N_vals(all_lost_lengths,total_removed_len);
            stat_output <<
                    "\nRejected Sequences (no frame detected):" <<
                    "\n\tTotal sequences: "                     << count_removed    <<
                    "\n\tAverage sequence length(bp): "         << avg_lost         <<
                    "\n\tn50: "                                 << removed_n.first  <<
                    "\n\tn90: "                                 << removed_n.second <<
                    "\n\tLongest sequence(bp): "  << max_removed<< " (" << max_removed_seq << ")" <<
                    "\n\tShortest sequence(bp): " << min_removed<< " (" << min_removed_seq << ")" <<"\n";
        }
        std::string stat_out_msg = stat_output.str();
        entapExecute::print_statistics(stat_out_msg,_outpath);
    } catch (ExceptionHandler &e) {throw e;}
    entapInit::print_msg("Success!");
}

FrameSelection::frame_map_type FrameSelection::genemark_parse_protein(std::string &protein) {
    frame_map_type protein_map;

    entapInit::print_msg("Parsing protein file at: " + protein);
    if (!entapInit::file_exists(protein)) {
        throw ExceptionHandler("File located at: " + protein + " does not exist!",
                               ENTAP_ERR::E_RUN_GENEMARK_PARSE);
    }
    std::string out_protein = protein + "alt";
    std::ofstream out_file(out_protein,std::ios::out | std::ios::app);
    std::ifstream in_file(protein);
    std::string line, sequence, seq_id;
    frame_seq protein_sequence;
    while (getline(in_file,line)){
        if (line.empty()) continue;
        if (line.find(">")==0) {
            if (!seq_id.empty()) {
                std::string sub = sequence.substr(sequence.find("\n")+1);
                long line_chars = std::count(sub.begin(),sub.end(),'\n');
                unsigned long seq_len = sub.length() - line_chars;
                protein_sequence = {seq_len,sequence, ""};
                protein_map.emplace(seq_id,protein_sequence);
            }
            unsigned long first = line.find(">")+1;
            unsigned long second = line.find("\t");
            seq_id = line.substr(first,second-first);
            sequence = line + "\n";
            out_file << ">" + seq_id << std::endl;
        } else {
            sequence += line + "\n";
            out_file << line << std::endl;
        }
    }
    out_file << line << std::endl;
    std::string sub = sequence.substr(sequence.find("\n")+1);
    long line_chars = std::count(sub.begin(),sub.end(),'\n');
    unsigned long seq_len = sub.length() - line_chars;
    protein_sequence = {seq_len,sequence, ""};
    protein_map.emplace(seq_id,protein_sequence);
    entapInit::print_msg("Success!");
    in_file.close(); out_file.close();
    boostFS::remove(protein);
    boostFS::rename(out_protein,protein);
    return protein_map;
}

void FrameSelection::genemark_parse_lst(std::string &lst_path, frame_map_type &current_map) {
    entapInit::print_msg("Parsing file at: " + lst_path);
    std::ifstream in_file(lst_path);
    std::string line, seq_id;
    while (getline(in_file,line)) {
        if (line.empty()) continue;
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        if (line.find("FASTA") == 0) {
            unsigned long first = line.find(":") + 1;
            seq_id = line.substr(first);
        } else if (isdigit(line.at(0))) {
            std::string frame;
            bool prime_5 = line.find("<") != std::string::npos;
            bool prime_3 = line.find(">") != std::string::npos;
            if (prime_5 && prime_3) {
                frame = "Internal";
            } else if (!prime_5 && !prime_3) {
                frame = "Complete";
            } else if (prime_5) {
                frame = ENTAP_EXECUTE::FRAME_SELECTION_FIVE_FLAG;
            } else frame = ENTAP_EXECUTE::FRAME_SELECTION_THREE_FLAG;
            std::map<std::string, frame_seq>::iterator it = current_map.find(seq_id);
            if (it != current_map.end()) {
                it->second.frame_type = frame;
            } else {
                throw ExceptionHandler("Sequence: " + seq_id + " not found in map",
                                       ENTAP_ERR::E_RUN_GENEMARK_PARSE);
            }
        }
    }
    entapInit::print_msg("Success!");
}



