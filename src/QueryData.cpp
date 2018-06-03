/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
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
#include "FileSystem.h"


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
QueryData::QueryData(std::string &input_file, std::string &out_path, UserInput *userinput,
                    FileSystem* filesystem) {
    FS_dprint("Processing transcriptome...");

    std::stringstream                        out_msg;
    std::string                              out_name;
    std::string                              out_new_path;
    std::string                              line;
    std::string                              sequence;
    std::string                              seq_id;
    std::string                              longest_seq;
    std::string                              shortest_seq;
    std::string                              transcript_type;
    uint32                                   count_seqs=0;
    uint64                                   total_len=0;
    uint16                                   shortest_len=10000;
    uint16                                   longest_len=0;
    uint16                                   len;
    fp64                                     avg_len;
    std::vector<uint16>                      sequence_lengths;
    std::pair<uint16, uint16>                n_vals;
    bool                                     is_complete;

    _total_sequences = 0;
    _pipeline_flags  = 0;
    _data_flags      = 0;
    _pSEQUENCES      = new QUERY_MAP_T;

    _pUserInput  = userinput;
    _pFileSystem = filesystem;

    _trim          = _pUserInput->has_input(UInput::INPUT_FLAG_TRIM);
    is_complete    = _pUserInput->has_input(UInput::INPUT_FLAG_COMPLETE);

    if (!_pFileSystem->file_exists(input_file)) {
        throw ExceptionHandler("Input transcriptome not found at: " + input_file,ERR_ENTAP_INPUT_PARSE);
    }

    out_name     = _pFileSystem->get_filename(input_file);
    out_new_path = PATHS(out_path,out_name);
    _pFileSystem->delete_file(out_new_path);

    set_input_type(input_file);
    _protein ? transcript_type = PROTEIN_FLAG : transcript_type = NUCLEO_FLAG;

    std::ifstream in_file(input_file);
    std::ofstream out_file(out_new_path,std::ios::out | std::ios::app);

    while (true) {
        std::getline(in_file, line);
        if (line.empty() && !in_file.eof()) continue;
        if (line.find(FASTA_FLAG) == 0 || in_file.eof()) {
            if (!seq_id.empty()) {
                if (in_file.eof()) {
                    out_file << line << std::endl;
                    sequence += line + "\n";
                }
                QuerySequence *query_seq = new QuerySequence(_protein,sequence, seq_id);
                if (is_complete) query_seq->setFrame(COMPLETE_FLAG);
                if (_pSEQUENCES->find(seq_id) != _pSEQUENCES->end()) {
                    throw ExceptionHandler("Duplicate headers in your input transcriptome: " + seq_id,
                        ERR_ENTAP_INPUT_PARSE);
                }
                _pSEQUENCES->emplace(seq_id, query_seq);
                count_seqs++;
                len = (uint16) query_seq->getSeq_length();
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
    _total_sequences = count_seqs;
    _protein  ? _start_prot_len = total_len : _start_nuc_len = total_len;
    // first - n50, second - n90
    n_vals = calculate_N_vals(sequence_lengths, total_len);

    out_msg <<std::fixed<<std::setprecision(2);
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
    _pFileSystem->print_stats(msg);
    FS_dprint("Success!");
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
        if (line.empty()) continue;
        if (line_count++ > LINE_COUNT) break;
        if (line_count < SEQ_DPRINT_CONUT) FS_dprint("Transcriptome Line: \n" + line);
        line.pop_back(); // Account for newline/other
        if (line.find('>') == std::string::npos) {
            for (char &c : line) {
                if (NUCLEO_MAP.find((char) toupper(c)) == NUCLEO_MAP.end()) deviations++;
            }
        }
    }
    _protein = deviations > NUCLEO_DEV;
    in_file.close();
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
std::pair<uint16, uint16> QueryData::calculate_N_vals
        (std::vector<uint16> &seq_lengths, uint64 total_len) {

    uint64 temp_len=0;
    uint16 n_50=0;
    uint16 n_90=0;
    fp64   fifty_len;
    fp64   ninety_len;

    std::sort(seq_lengths.begin(),seq_lengths.end(), std::greater<uint16>());
    fifty_len  = total_len * N_50_PERCENT;
    ninety_len = total_len * N_90_PERCENT;
    for (uint16 val : seq_lengths) {
        temp_len += val;
        if (temp_len > fifty_len && n_50 == 0) n_50 = val;
        if (temp_len > ninety_len) {
            n_90 = val;
            break;
        }
    }
    return std::pair<uint16, uint16> (n_50,n_90);
}

// Will be used later
std::pair<uint16, uint16> QueryData::calculate_N_vals (ExecuteStates state, bool kept) {
    std::vector<uint16> seq_len_vect;
    uint64 total_len=0;     // Kept sequences only
    uint32 total_seq=0;     // Kept sequences only
    uint64 temp_len=0;
    uint64 n_50=0;
    uint64 n_90=0;
    fp64   fifty_len;
    fp64   ninety_len;

    // Recalculate based upon what sequences are left (kept)
    for (auto &pair : *_pSEQUENCES) {
        if (pair.second->is_kept()) {
            total_len += pair.second->getSeq_length();
            total_seq++;
            seq_len_vect.push_back((uint16)pair.second->getSeq_length());
        }
    }

    std::sort(seq_len_vect.begin(),seq_len_vect.end());
    fifty_len  = total_len * N_50_PERCENT;
    ninety_len = total_len * N_90_PERCENT;

    for (uint16 &val : seq_len_vect) {
        temp_len += val;
        if (temp_len > fifty_len && n_50 == 0) n_50 = val;
        if (temp_len > ninety_len) {
            n_90 = val;
            break;
        }
    }

    // Print new transcriptome stats

    return std::pair<uint16, uint16> (n_50,n_90);
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
    for (auto &pair : *_pSEQUENCES) {
        switch (state) {
            case EXPRESSION_FILTERING:
                pair.second->QUERY_FLAG_SET(QuerySequence::QUERY_EXPRESSION_KEPT);
                break;
            case FRAME_SELECTION:
                pair.second->QUERY_FLAG_SET(QuerySequence::QUERY_IS_PROTEIN);        // Probably already done
                pair.second->QUERY_FLAG_SET(QuerySequence::QUERY_FRAME_KEPT);
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
void QueryData::final_statistics(std::string &outpath, std::vector<uint16> &ontology_flags) {
    FS_dprint("Pipeline finished! Calculating final statistics...");

    std::stringstream      ss;
    uint32                 count_total_sequences=0;
    uint32                 count_exp_kept=0;
    uint32                 count_exp_reject=0;
    uint32                 count_frame_kept=0;
    uint32                 count_frame_rejected=0;
    uint32                 count_sim_hits=0;
    uint32                 count_sim_no_hits=0;
    uint32                 count_ontology=0;
    uint32                 count_no_ontology=0;
    uint32                 count_one_go=0;
    uint32                 count_one_kegg=0;
    uint32                 count_sim_only=0;
    uint32                 count_ontology_only=0;
    uint32                 count_TOTAL_ann=0;
    uint32                 count_TOTAL_unann=0;
    std::string            out_unannotated_nucl_path;
    std::string            out_unannotated_prot_path;
    std::string            out_annotated_nucl_path;
    std::string            out_annotated_prot_path;
    std::string            out_msg;
    bool                   is_exp_kept;
    bool                   is_prot;
    bool                   is_hit;
    bool                   is_ontology;
    bool                   is_one_go;
    bool                   is_one_kegg;

    out_unannotated_nucl_path = PATHS(outpath, OUT_UNANNOTATED_NUCL);
    out_unannotated_prot_path = PATHS(outpath, OUT_UNANNOTATED_PROT);
    out_annotated_nucl_path   = PATHS(outpath, OUT_ANNOTATED_NUCL);
    out_annotated_prot_path   = PATHS(outpath, OUT_ANNOTATED_PROT);

    // Re-write these files
    _pFileSystem->delete_file(out_unannotated_nucl_path);
    _pFileSystem->delete_file(out_unannotated_prot_path);
    _pFileSystem->delete_file(out_annotated_nucl_path);
    _pFileSystem->delete_file(out_annotated_prot_path);

    std::ofstream file_unannotated_nucl(out_unannotated_nucl_path, std::ios::out | std::ios::app);
    std::ofstream file_unannotated_prot(out_unannotated_prot_path, std::ios::out | std::ios::app);
    std::ofstream file_annotated_nucl(out_annotated_nucl_path, std::ios::out | std::ios::app);
    std::ofstream file_annotated_prot(out_annotated_prot_path, std::ios::out | std::ios::app);

    for (auto &pair : *_pSEQUENCES) {
        count_total_sequences++;
        is_exp_kept = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_EXPRESSION_KEPT);
        is_prot = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN);
        is_hit = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_BLAST_HIT);
        is_ontology = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_FAMILY_ASSIGNED); // TODO Fix for interpro
        is_one_go = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_ONE_GO);
        is_one_kegg = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_ONE_KEGG);

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
            if (!pair.second->get_sequence_n().empty())
                file_annotated_nucl<<pair.second->get_sequence_n()<<std::endl;
            if (!pair.second->get_sequence_p().empty()) {
                file_annotated_prot<<pair.second->get_sequence_p()<<std::endl;
            }
        } else {
            // Not annotated
            if (!pair.second->get_sequence_n().empty())
                file_unannotated_nucl<<pair.second->get_sequence_n()<<std::endl;
            if (!pair.second->get_sequence_p().empty()) {
                file_unannotated_prot<<pair.second->get_sequence_p()<<std::endl;
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

    if (DATA_FLAG_GET(SUCCESS_EXPRESSION)) {
        ss <<
           "\nExpression Analysis" <<
           "\n\tKept sequences: "  << count_exp_kept    <<
           "\n\tLost sequences: "  << count_exp_reject;
    }
    if (DATA_FLAG_GET(SUCCESS_FRAME_SEL)) {
        ss <<
           "\nFrame Selection"              <<
           "\n\tTotal sequences retained: " << count_frame_kept     <<
           "\n\tTotal sequences removed: "  << count_frame_rejected;
    }
    if (DATA_FLAG_GET(SUCCESS_SIM_SEARCH)) {
        ss <<
           "\nSimilarity Search"                               <<
           "\n\tTotal unique sequences with an alignment: "    << count_sim_hits <<
           "\n\tTotal unique sequences without an alignment: " << count_sim_no_hits;
    }
    if (DATA_FLAG_GET(SUCCESS_ONTOLOGY)) {
        for (uint16 flag : ontology_flags) {
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
    _pFileSystem->print_stats(out_msg);
}

std::string QueryData::trim_sequence_header(std::string &header, std::string line) {
    std::string   sequence;
    int16         pos;

    if (line.find(">") != std::string::npos) {
        pos = (int16) line.find(">");
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

QUERY_MAP_T* QueryData::get_sequences_ptr() {
    return this->_pSEQUENCES;
}

QueryData::~QueryData() {
    FS_dprint("Killing object - QueryData");
    for(QUERY_MAP_T::iterator it = _pSEQUENCES->begin(); it != _pSEQUENCES->end(); it++) {
        delete it->second;
        it->second = nullptr;
    }
    FS_dprint("QuerySequence data freed");
    SAFE_DELETE(_pSEQUENCES);
}

bool QueryData::DATA_FLAG_GET(DATA_FLAGS flag) {
    return (_data_flags & flag) != 0;
}

void QueryData::DATA_FLAG_SET(DATA_FLAGS flag) {
    _data_flags |= flag;
}

void QueryData::DATA_FLAG_CLEAR(DATA_FLAGS flag) {
    _data_flags &= ~flag;
}

QuerySequence *QueryData::get_sequence(std::string &query_id) {
    QUERY_MAP_T::iterator it = _pSEQUENCES->find(query_id);
    if (it != _pSEQUENCES->end()) {
        // Sequence found, return
        return it->second;
    } else {
        // Sequence NOT found retun null
        return nullptr;
    }
}


