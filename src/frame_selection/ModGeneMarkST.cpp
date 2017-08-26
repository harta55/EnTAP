/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/


#include <iomanip>
#include "ModGeneMarkST.h"
#include "../ExceptionHandler.h"
#include "../EntapExecute.h"

std::pair<bool, std::string> ModGeneMarkST::verify_files() {
    std::string lst_file;

    boost::filesystem::path file_name(boostFS::path(_inpath).filename());
    _final_out = (boostFS::path(_frame_outpath) / file_name).string() + ".faa";
    lst_file = file_name.string() + ".lst";
    _final_lst = (boostFS::path(_frame_outpath) / boostFS::path(lst_file)).string();
    if (file_exists(_final_out) && file_exists(_final_lst)) {
        print_debug("File found at: " + _final_out + "\n"
                "continuing enTAP with this file and skipping frame selection");
        return std::make_pair(true, _final_out);
    }
    print_debug("File not found at " + _final_out + " so continuing frame selection");
    return std::make_pair(false, "");
}

std::string ModGeneMarkST::execute(std::map<std::string, QuerySequence> &SEQUENCES) {
    // Outfiles: file/path.faa, file/path.fnn
    // assumes working directory as output right now
    std::string     lst_file;
    std::string     out_gmst_log;
    std::string     out_hmm_file;
    std::string     genemark_cmd;
    std::string     genemark_std_out;
    std::string     line;

    boost::filesystem::path file_name(boostFS::path(_inpath).filename());
    std::list<std::string> out_names {file_name.string()+".faa",
                                      file_name.string()+".fnn"};
    lst_file = file_name.string() + ".lst";
    out_gmst_log = (boostFS::path(_frame_outpath) / boostFS::path(GENEMARK_LOG_FILE)).string();
    out_hmm_file = (boostFS::path(_frame_outpath) / boostFS::path(GENEMARK_HMM_FILE)).string();

    genemark_cmd = _exe_path + " -faa -fnn " + _inpath;
    genemark_std_out = (boostFS::path(_frame_outpath) / boostFS::path(GENEMARK_STD_OUT)).string();
    print_debug("Running genemark...\n" + genemark_cmd);

    if (execute_cmd(genemark_cmd,genemark_std_out) != 0 ) {
        throw ExceptionHandler("Error in running genemark at file located at: " +
                               _inpath, ENTAP_ERR::E_INIT_INDX_DATA_NOT_FOUND);
    }
    print_debug("Success!");
    // Format genemarks-t output (remove blank lines)
    print_debug("Formatting genemark files");

    for (std::string path : out_names) {
        std::ifstream in_file(path);
        std::string temp_name = path+"_alt";
        std::string out_path = (boostFS::path(_frame_outpath) / boostFS::path(path)).string();
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
    if (rename(lst_file.c_str(),_final_lst.c_str())!=0 ||
        rename(GENEMARK_LOG_FILE.c_str(),out_gmst_log.c_str())!=0 ) {
        throw ExceptionHandler("Error moving genemark results", ENTAP_ERR::E_INIT_TAX_READ);
    }
    if (file_exists(GENEMARK_HMM_FILE)) {
        rename(GENEMARK_HMM_FILE.c_str(),out_hmm_file.c_str());
    }
    print_debug("Success!");
    return _final_out;
}


/**
 *
 * @param protein_path - path to .faa file
 * @param lst_path - path to .lst genemark file
 */
void ModGeneMarkST::parse(std::map<std::string, QuerySequence> &SEQUENCES) {
    // generate maps, query->sequence
    print_debug("Beginning to calculate Genemark statistics...");

    std::string                             out_removed_path;
    std::string                             out_internal_path;
    std::string                             out_complete_path;
    std::string                             out_partial_path;
    std::string                             figure_results_path;
    std::string                             figure_results_png;
    std::string                             figure_removed_path;
    std::string                             figure_removed_png;
    std::string                             min_removed_seq;
    std::string                             min_kept_seq;
    std::string                             max_removed_seq;
    std::string                             max_kept_seq;
    std::stringstream                       stat_output;
    std::map<std::string, std::ofstream*>   file_map;
    std::map<std::string, unsigned long>    count_map;
    std::vector<unsigned long>              all_kept_lengths;
    std::vector<unsigned long>              all_lost_lengths;
    float                                   avg_selected;
    float                                   avg_lost;
    std::pair<unsigned long, unsigned long> kept_n;
    GraphingStruct                          graphingStruct;

    boostFS::remove_all(_processed_path);
    boostFS::remove_all(_figure_path);
    boostFS::create_directories(_processed_path);
    boostFS::create_directories(_figure_path);
    out_removed_path    = (boostFS::path(_processed_path) / boostFS::path(FRAME_SELECTION_LOST)).string();
    out_internal_path   = (boostFS::path(_processed_path) / boostFS::path(FRAME_SELECTION_INTERNAL)).string();
    out_complete_path   = (boostFS::path(_processed_path) / boostFS::path(FRAME_SELECTION_COMPLTE)).string();
    out_partial_path    = (boostFS::path(_processed_path) / boostFS::path(FRAME_SELECTION_PARTIAL)).string();
    figure_removed_path = (boostFS::path(_figure_path)    / boostFS::path(GRAPH_TEXT_REF_COMPAR)).string();
    figure_removed_png  = (boostFS::path(_figure_path)    / boostFS::path(GRAPH_FILE_REF_COMPAR)).string();
    figure_results_path = (boostFS::path(_figure_path)    / boostFS::path(GRAPH_TEXT_FRAME_RESUTS)).string();
    figure_results_png  = (boostFS::path(_figure_path)    / boostFS::path(GRAPH_FILE_FRAME_RESUTS)).string();

    // all nucleotide lengths
    unsigned long min_removed=10000,min_selected=10000,max_removed=0,max_selected=0,
            total_removed_len=0,total_kept_len=0, count_selected = 0, count_removed = 0,
            count_partial_5=0, count_partial_3=0, count_internal=0, count_complete=0;
    try {
        std::map<std::string,frame_seq> protein_map = genemark_parse_protein(_final_out);
        genemark_parse_lst(_final_lst,protein_map);

        std::ofstream file_figure_removed(figure_removed_path,std::ios::out | std::ios::app);
        std::ofstream file_figure_results(figure_results_path,std::ios::out | std::ios::app);

        file_figure_removed << "flag\tsequence length" << std::endl;    // First line placeholder, not used
        file_figure_results << "flag\tsequence length" << std::endl;

        file_map[FRAME_SELECTION_LOST_FLAG] =
                new std::ofstream(out_removed_path, std::ios::out | std::ios::app);
        file_map[FRAME_SELECTION_INTERNAL_FLAG] =
                new std::ofstream(out_internal_path, std::ios::out | std::ios::app);
        file_map[FRAME_SELECTION_COMPLETE_FLAG] =
                new std::ofstream(out_complete_path, std::ios::out | std::ios::app);
        file_map[FRAME_SELECTION_FIVE_FLAG] =
                new std::ofstream(out_partial_path, std::ios::out | std::ios::app);
        file_map[FRAME_SELECTION_THREE_FLAG] = file_map[FRAME_SELECTION_FIVE_FLAG];

        count_map ={
                {FRAME_SELECTION_INTERNAL_FLAG,count_internal},
                {FRAME_SELECTION_COMPLETE_FLAG,count_complete},
                {FRAME_SELECTION_FIVE_FLAG,count_partial_5},
                {FRAME_SELECTION_THREE_FLAG,count_partial_3},
        };

        for (auto& pair : SEQUENCES) {
            std::map<std::string,frame_seq>::iterator p_it = protein_map.find(pair.first);
            if (!pair.second.is_is_expression_kept()) continue; // Skip seqs that were lost to expression
            if (p_it != protein_map.end()) {
                // Kept sequence, either partial, complete, or internal
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
                file_figure_removed << GRAPH_KEPT_FLAG << '\t' << std::to_string(length) << std::endl;
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
                *file_map[FRAME_SELECTION_LOST_FLAG] << pair.second.get_sequence_n() << std::endl;
                unsigned long length = pair.second.getSeq_length();  // Nucleotide sequence length

                if (length < min_removed) {
                    min_removed = length;
                    min_removed_seq = pair.first;
                }
                if (length > max_removed) {
                    max_removed_seq = pair.first;
                    max_removed = length;
                }
                file_figure_removed << GRAPH_REJECTED_FLAG << '\t' << std::to_string(length) << std::endl;
                all_lost_lengths.push_back(length);
                total_removed_len += length;
            }
        }
        for(auto& pair : file_map) {
            if (pair.first.compare(FRAME_SELECTION_THREE_FLAG)!=0){
                pair.second->close();
                delete pair.second;
                pair.second = 0;
            }
        }
        // Calculate and print stats
        avg_selected = (float)total_kept_len / count_selected;
        avg_lost = (float)total_removed_len / count_removed;
        stat_output<<std::fixed<<std::setprecision(2);
        stat_output <<
                    ENTAP_STATS::SOFTWARE_BREAK             <<
                    "Frame Selection: GenemarkS-T\n"        <<
                    ENTAP_STATS::SOFTWARE_BREAK             <<
                    "Total sequences frame selected: "      << count_selected          <<
                    "\n\tTranslated protein sequences: "    << _final_out              <<
                    "\nTotal sequences removed (no frame): "<< count_removed           <<
                    "\n\tFrame selected CDS removed: "      << out_removed_path        <<
                    "\nTotal of "                           <<
                    count_map[FRAME_SELECTION_FIVE_FLAG]    << " 5 prime partials and "<<
                    count_map[FRAME_SELECTION_THREE_FLAG]   << " 3 prime partials"     <<
                    "\n\tPartial CDS: "                     << out_partial_path        <<
                    "\nTotal of "                           <<
                    count_map[FRAME_SELECTION_COMPLETE_FLAG]<<" complete genes:\n\t" << out_complete_path<<
                    "\nTotal of "                           <<
                    count_map[FRAME_SELECTION_INTERNAL_FLAG]<<" internal genes:\n\t" << out_internal_path<<"\n\n";

        stat_output <<
                    ENTAP_STATS::SOFTWARE_BREAK                                <<
                    "Frame Selection: New Reference Transcriptome Statistics\n"<<
                    ENTAP_STATS::SOFTWARE_BREAK;

        kept_n = entapExecute::calculate_N_vals(all_kept_lengths,total_kept_len);
        stat_output <<
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
                        "\nRemoved Sequences (no frame):"           <<
                        "\n\tTotal sequences: "                     << count_removed    <<
                        "\n\tAverage sequence length(bp): "         << avg_lost         <<
                        "\n\tn50: "                                 << removed_n.first  <<
                        "\n\tn90: "                                 << removed_n.second <<
                        "\n\tLongest sequence(bp): "  << max_removed<< " (" << max_removed_seq << ")" <<
                        "\n\tShortest sequence(bp): " << min_removed<< " (" << min_removed_seq << ")" <<"\n";
        }
        std::string stat_out_msg = stat_output.str();
        print_statistics(stat_out_msg);

        // Figure handling
        file_figure_results << GRAPH_REJECTED_FLAG           << '\t' << std::to_string(count_removed)   <<std::endl;
        file_figure_results << FRAME_SELECTION_FIVE_FLAG     << '\t' << std::to_string(count_map[FRAME_SELECTION_FIVE_FLAG]) <<std::endl;
        file_figure_results << FRAME_SELECTION_THREE_FLAG    << '\t' << std::to_string(count_map[FRAME_SELECTION_THREE_FLAG]) <<std::endl;
        file_figure_results << FRAME_SELECTION_COMPLETE_FLAG << '\t' << std::to_string(count_map[FRAME_SELECTION_COMPLETE_FLAG])   <<std::endl;
        file_figure_results << FRAME_SELECTION_INTERNAL_FLAG << '\t' << std::to_string(count_map[FRAME_SELECTION_INTERNAL_FLAG])   <<std::endl;

        graphingStruct.text_file_path = figure_results_path;
        graphingStruct.graph_title    = GRAPH_TITLE_FRAME_RESULTS;
        graphingStruct.fig_out_path   = figure_results_png;
        graphingStruct.software_flag  = GRAPH_FRAME_FLAG;
        graphingStruct.graph_type     = GRAPH_PIE_RESULTS_FLAG;
        pGraphingManager->graph(graphingStruct);

        graphingStruct.text_file_path = figure_removed_path;
        graphingStruct.graph_title    = GRAPH_TITLE_REF_COMPAR;
        graphingStruct.fig_out_path   = figure_removed_png;
        graphingStruct.graph_type     = GRAPH_COMP_BOX_FLAG;
        pGraphingManager->graph(graphingStruct);

    } catch (ExceptionHandler &e) {throw e;}
    print_debug("Success!");
}


ModGeneMarkST::frame_map_t ModGeneMarkST::genemark_parse_protein(std::string &protein) {
    print_debug("Parsing protein file at: " + protein);

    frame_map_t     protein_map;
    std::string     out_protein;
    frame_seq       protein_sequence;
    std::string     line;
    std::string     sequence;
    std::string     seq_id;

    if (!file_exists(protein)) {
        throw ExceptionHandler("File located at: " + protein + " does not exist!",
                               ENTAP_ERR::E_RUN_GENEMARK_PARSE);
    }
    out_protein = protein + "alt";
    std::ofstream out_file(out_protein,std::ios::out | std::ios::app);
    std::ifstream in_file(protein);
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
    print_debug("Success!");
    in_file.close(); out_file.close();
    boostFS::remove(protein);
    boostFS::rename(out_protein,protein);
    return protein_map;
}

void ModGeneMarkST::genemark_parse_lst(std::string &lst_path, frame_map_t &current_map) {
    print_debug("Parsing file at: " + lst_path);

    std::string     line;
    std::string     seq_id;
    bool            prime_5;
    bool            prime_3;

    std::ifstream in_file(lst_path);
    while (getline(in_file,line)) {
        if (line.empty()) continue;
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        if (line.find("FASTA") == 0) {
            unsigned long first = line.find(":") + 1;
            seq_id = line.substr(first);
        } else if (isdigit(line.at(0))) {
            std::string frame;
            prime_5 = line.find("<") != std::string::npos;
            prime_3 = line.find(">") != std::string::npos;
            if (prime_5 && prime_3) {
                frame = "Internal";
            } else if (!prime_5 && !prime_3) {
                frame = "Complete";
            } else if (prime_5) {
                frame = FRAME_SELECTION_FIVE_FLAG;
            } else frame = FRAME_SELECTION_THREE_FLAG;
            std::map<std::string, frame_seq>::iterator it = current_map.find(seq_id);
            if (it != current_map.end()) {
                it->second.frame_type = frame;
            } else {
                throw ExceptionHandler("Sequence: " + seq_id + " not found in map",
                                       ENTAP_ERR::E_RUN_GENEMARK_PARSE);
            }
        }
    }
    print_debug("Success!");
}