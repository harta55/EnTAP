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

#include <iomanip>
#include "QueryData.h"
#include "ExceptionHandler.h"


/**
 * ======================================================================
 * Function QueryData::QueryData(std::string &input_file,
 *                              std::string &out_path,
 *                              bool &is_complete,
 *                              bool &trim)
 *
 * Description          - Parses input transcriptome and converts to map of
 *                        each query sequence
 *                      - This map is passed throughout EnTAP execution and
 *                        updated
 *
 * Notes                - None
 *
 * @param input_file    - Path to input transcriptome
 * @param trim          - Flag from user to trim sequence ID to first space
 * @param is_complete   - Flag from user if the entire transcriptome is a
 *                        complete gene
 * @return              - None
 *
 * =====================================================================
 */
QueryData::QueryData(std::string &input_file, std::string &out_path, bool &is_complete, bool &trim) {
    print_debug("Processing transcriptome...");

    std::stringstream                        out_msg;
    std::string                              out_name;
    std::string                              out_new_path;
    std::string                              line;
    std::string                              sequence;
    std::string                              seq_id;
    std::string                              longest_seq;
    std::string                              shortest_seq;
    std::string                              transcript_type;
    unsigned long                            count_seqs=0;
    unsigned long long                       total_len=0;
    unsigned short                           shortest_len=10000;
    unsigned short                           longest_len=0;
    unsigned short                           len;
    double                                   avg_len;
    std::vector<unsigned short>              sequence_lengths;
    std::pair<unsigned short, unsigned short>n_vals;

    _trim = trim;

    if (!file_exists(input_file)) {
        throw ExceptionHandler("Input file not found at: " + input_file,ENTAP_ERR::E_INPUT_PARSE);
    }

    boostFS::path path(input_file);
    out_name     = path.filename().string();
    out_new_path = PATHS(out_path,out_name);
    boostFS::remove(out_new_path);

    set_input_type(input_file);
    _protein ? transcript_type = PROTEIN_FLAG : transcript_type = NUCLEO_FLAG;

    std::ifstream in_file(input_file);
    std::ofstream out_file(out_new_path,std::ios::out | std::ios::app);

    while (true) {
        std::getline(in_file, line);
        if (line.empty() && !in_file.eof()) continue;
        if (line.find(">") == 0 || in_file.eof()) {
            if (!seq_id.empty()) {
                if (in_file.eof()) {
                    out_file << line << std::endl;
                    sequence += line + "\n";
                }
                QuerySequence query_seq = QuerySequence(_protein,sequence);
                if (is_complete) query_seq.setFrame(COMPLETE_FLAG);
                query_seq.setQseqid(seq_id);
                _SEQUENCES.emplace(seq_id, query_seq);
                count_seqs++;
                len = (unsigned short) query_seq.getSeq_length();
                total_len += len;
                if (len > longest_len) {
                    longest_len = len;longest_seq = seq_id;
                }
                if (len < shortest_len) {
                    shortest_len = len;shortest_seq = seq_id;
                }
                sequence_lengths.push_back(len);
            }
            if (in_file.eof()) break;
            sequence = trim_sequence_header(seq_id, line);
            out_file << sequence;
        } else {
            out_file << line << std::endl;
            sequence += line + "\n";
        }
    }
    in_file.close();
    out_file.close();
    avg_len = total_len / count_seqs;
    // first - n50, second - n90
    n_vals = calculate_N_vals(sequence_lengths, total_len);

    out_msg<<std::fixed<<std::setprecision(2);
    out_msg << ENTAP_STATS::SOFTWARE_BREAK
            << "Transcriptome Statistics\n"
            << ENTAP_STATS::SOFTWARE_BREAK<<
            transcript_type << " sequences found"          <<
            "\nTotal sequences: "                          << count_seqs    <<
            "\nTotal length of transcriptome(bp): "        << total_len     <<
            "\nAverage sequence length(bp): "              << avg_len       <<
            "\nn50: "                                      << n_vals.first  <<
            "\nn90: "                                      << n_vals.second <<
            "\nLongest sequence(bp): " << longest_len << " ("<<longest_seq<<")"<<
            "\nShortest sequence(bp): "<< shortest_len<<" ("<<shortest_seq<<")";
    if (is_complete)out_msg<<"\nAll sequences ("<<count_seqs<<") were flagged as complete genes";
    std::string msg = out_msg.str();
    print_statistics(msg);
    print_debug("Success!");
    input_file = out_new_path;
}


void QueryData::set_input_type(std::string &in) {
    std::string    line;
    uint8          line_count;
    uint16         deviations;
    std::ifstream in_file(in);

    line_count = 0;
    deviations = 0;
    while(std::getline(in_file,line)) {
        if (line.empty() || line_count > LINE_COUNT) continue;
        if (line.find('>') == std::string::npos) {
            for (char &c : line) {
                toupper(c);
                if (NUCLEO_MAP.find(c) == NUCLEO_MAP.end()) deviations++;
            }
        }
        line_count++;
    }
    _protein = deviations > NUCLEO_DEV;
    in_file.close();
}

bool QueryData::is_protein() const {
    return _protein;
}


/**
 * Description - This function calculates n50 and n90 values with sequence
 *               length (nucleotide) information from a transcriptome
 *
 * @param seq_lengths - Vector of all (nucl) sequence lengths in transcriptome.
 *                    These will be sorted.
 * @param total_len   - Sum of all nucleotide lengths in transcriptome
 *
 * @return            - Pair of <n50,n90>
 */
std::pair<unsigned short, unsigned short> QueryData::calculate_N_vals
        (std::vector<unsigned short> &seq_lengths, unsigned long long total_len) {
    std::sort(seq_lengths.begin(),seq_lengths.end());
    unsigned long temp_len=0, n_50=0,n_90=0;
    double fifty_len = total_len * 0.5;
    double ninety_len = total_len * 0.9;
    for (unsigned long val : seq_lengths) {
        temp_len += val;
        if (temp_len > fifty_len && n_50 == 0) n_50 = val;
        if (temp_len > ninety_len) {
            n_90 = val;
            break;
        }
    }
    return std::pair<unsigned short, unsigned short> (n_50,n_90);
}


/**
 * ======================================================================
 * Function void flag_transcripts(ExecuteStates state,
 *                              std::map<std::string, QuerySequence>& map)
 *
 * Description          - Sets boolean flags if a certain stage is skipped
 *                        in pipeline
 *                      - Used for statistics
 *
 * Notes                - Will be moved to transcriptome/project object
 *
 * @param state         - Current enum state of execution
 *
 * @return              - None
 * ======================================================================
 */
void QueryData::flag_transcripts(ExecuteStates state) {
    for (auto &pair : _SEQUENCES) {
        switch (state) {
            case RSEM:
                pair.second.set_is_expression_kept(true);
                break;
            case FRAME_SELECTION:
                pair.second.setIs_protein(true);        // Probably already done
                break;
            default:
                break;
        }
    }
}



/**
 * ======================================================================
 * Function final_statistics(std::map<std::string, QuerySequence> &SEQUENCE_MAP)
 *
 * Description          - Calculates final statistical information after
 *                        completed execution
 *                      - Compiles stats on each stage of pipeline
 *
 * Notes                - None
 *
 * @param SEQUENCE_MAP  - Map of each query sequence + data
 *
 * @return              - None
 *
 * =====================================================================
 */
void QueryData::final_statistics(std::string &outpath, std::vector<uint8> &ontology_flags) {
    print_debug("Pipeline finished! Calculating final statistics...");

    std::stringstream      ss;
    unsigned int           count_total_sequences=0;
    unsigned int           count_exp_kept=0;
    unsigned int           count_exp_reject=0;
    unsigned int           count_frame_kept=0;
    unsigned int           count_frame_rejected=0;
    unsigned int           count_sim_hits=0;
    unsigned int           count_sim_no_hits=0;
    unsigned int           count_ontology=0;
    unsigned int           count_no_ontology=0;
    unsigned int           count_one_go=0;
    unsigned int           count_one_kegg=0;
    unsigned int           count_sim_only=0;
    unsigned int           count_ontology_only=0;
    unsigned int           count_TOTAL_ann=0;
    unsigned int           count_TOTAL_unann=0;
    std::string            out_unannotated_nucl_path;
    std::string            out_unannotated_prot_path;
    std::string            out_annotated_nucl_path;
    std::string            out_annotated_prot_path;
    std::string            out_msg;

    out_unannotated_nucl_path = PATHS(outpath, OUT_UNANNOTATED_NUCL);
    out_unannotated_prot_path = PATHS(outpath, OUT_UNANNOTATED_PROT);
    out_annotated_nucl_path   = PATHS(outpath, OUT_ANNOTATED_NUCL);
    out_annotated_prot_path   = PATHS(outpath, OUT_ANNOTATED_PROT);

    // Switch to entap_out dir
    boostFS::remove(out_unannotated_nucl_path);
    boostFS::remove(out_unannotated_prot_path);
    boostFS::remove(out_annotated_nucl_path);
    boostFS::remove(out_annotated_prot_path);

    std::ofstream file_unannotated_nucl(out_unannotated_nucl_path, std::ios::out | std::ios::app);
    std::ofstream file_unannotated_prot(out_unannotated_prot_path, std::ios::out | std::ios::app);
    std::ofstream file_annotated_nucl(out_annotated_nucl_path, std::ios::out | std::ios::app);
    std::ofstream file_annotated_prot(out_annotated_prot_path, std::ios::out | std::ios::app);

    for (auto &pair : _SEQUENCES) {
        count_total_sequences++;
        bool is_exp_kept = pair.second.is_is_expression_kept();
        bool is_prot = pair.second.isIs_protein();
        bool is_hit = pair.second.is_is_database_hit();
        bool is_ontology = pair.second.is_is_family_assigned(); // TODO Fix for interpro
        bool is_one_go = pair.second.is_is_one_go();
        bool is_one_kegg = pair.second.is_is_one_kegg();

        is_exp_kept ? count_exp_kept++ : count_exp_reject++;
        is_prot ? count_frame_kept++ : count_frame_rejected++;
        is_hit ? count_sim_hits++ : count_sim_no_hits++;
        is_ontology ? count_ontology++ : count_no_ontology++;
        if (is_one_go) count_one_go++;
        if (is_one_kegg) count_one_kegg++;

        if (is_hit && !is_ontology) count_sim_only++;
        if (!is_hit && is_ontology) count_ontology_only++;

        if (is_hit || is_ontology) {
            // Is annotated
            count_TOTAL_ann++;
            if (!pair.second.get_sequence_n().empty())
                file_annotated_nucl<<pair.second.get_sequence_n()<<std::endl;
            if (!pair.second.get_sequence_p().empty()) {
                file_annotated_prot<<pair.second.get_sequence_p()<<std::endl;
            }
        } else {
            // Not annotated
            if (!pair.second.get_sequence_n().empty())
                file_unannotated_nucl<<pair.second.get_sequence_n()<<std::endl;
            if (!pair.second.get_sequence_p().empty()) {
                file_unannotated_prot<<pair.second.get_sequence_p()<<std::endl;
            }
            count_TOTAL_unann++;
        }
    }

    file_unannotated_nucl.close();
    file_unannotated_prot.close();
    file_annotated_nucl.close();
    file_annotated_prot.close();

    ss <<
       ENTAP_STATS::SOFTWARE_BREAK          <<
       "Final Annotation Statistics\n"      <<
       ENTAP_STATS::SOFTWARE_BREAK          <<
       "Total Sequences: "                  << count_total_sequences;

    if (_EXPRESSION_SUCCESS) {
        ss <<
           "\nExpression Analysis" <<
           "\n\tKept sequences: "  << count_exp_kept    <<
           "\n\tLost sequences: "  << count_exp_reject;
    }
    if (_FRAME_SELECTION_SUCCESS) {
        ss <<
           "\nFrame Selection"              <<
           "\n\tTotal sequences retained: " << count_frame_kept     <<
           "\n\tTotal sequences removed: "  << count_frame_rejected;
    }
    if (_SIM_SEARCH_SUCCESS) {
        ss <<
           "\nSimilarity Search"                               <<
           "\n\tTotal unique sequences with an alignment: "    << count_sim_hits <<
           "\n\tTotal unique sequences without an alignment: " << count_sim_no_hits;
    }
    if (_ONTOLOGY_SUCCESS) {
        for (uint8 flag : ontology_flags) {
            switch (flag) {
                case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
                    ss <<
                       "\nGene Families"        <<
                       "\n\tTotal unique sequences with family assignment: "    << count_ontology   <<
                       "\n\tTotal unique sequences without family assignment: " << count_no_ontology<<
                       "\n\tTotal unique sequences with at least one GO term: " << count_one_go     <<
                       "\n\tTotal unique sequences with at least one pathway (KEGG) assignment: "   << count_one_kegg;
                    break;
                case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
                    ss <<
                       "\nFinal InterPro stats coming soon!";
                    break;
                default:
                    break;
            }
        }
    }
    ss <<
       "\nTotals"   <<
       "\n\tTotal unique sequences annotated (similarity search alignments only): "      << count_sim_only      <<
       "\n\tTotal unique sequences annotated (gene family assignment only): "            << count_ontology_only <<
       "\n\tTotal unique sequences annotated (gene family and/or similarity search): "   << count_TOTAL_ann     <<
       "\n\tTotal unique sequences unannotated (gene family and/or similarity search): " << count_TOTAL_unann;

    out_msg = ss.str();
    print_statistics(out_msg);
}

std::string QueryData::trim_sequence_header(std::string &header, std::string line) {
    std::string   sequence;
    signed short  pos;

    if (line.find(">") != std::string::npos) {
        pos = (signed short) line.find(">");
    } else pos = -1;
    if (_trim) {
        if (line.find(" ") != std::string::npos) {
            header = line.substr(pos+1, line.find(" ")-1);
        } else header = line.substr(pos+1);
        sequence = ">" + header + "\n";
    } else {
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        header = line.substr(pos+1);
        sequence = line + "\n";
    }
    return sequence;
}

void QueryData::set_frame_stats(const FrameStats &_frame_stats) {
    QueryData::_frame_stats = _frame_stats;
}

void QueryData::set_EXPRESSION_SUCCESS(bool _EXPRESSION_SUCCESS) {
    QueryData::_EXPRESSION_SUCCESS = _EXPRESSION_SUCCESS;
}

void QueryData::set_FRAME_SELECTION_SUCCESS(bool _FRAME_SELECTION_SUCCESS) {
    QueryData::_FRAME_SELECTION_SUCCESS = _FRAME_SELECTION_SUCCESS;
}

void QueryData::set_SIM_SEARCH_SUCCESS(bool _SIM_SEARCH_SUCCESS) {
    QueryData::_SIM_SEARCH_SUCCESS = _SIM_SEARCH_SUCCESS;
}

void QueryData::set_ONTOLOGY_SUCCESS(bool _ONTOLOGY_SUCCESS) {
    QueryData::_ONTOLOGY_SUCCESS = _ONTOLOGY_SUCCESS;
}

QUERY_MAP_T* QueryData::get_sequences_ptr() {
    return &this->_SEQUENCES;
}
