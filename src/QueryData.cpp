/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2020, Alexander Hart, Dr. Jill Wegrzyn
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

#include "QueryData.h"
#include "QueryAlignment.h"
#include "ExceptionHandler.h"
#include "FileSystem.h"
#include "UserInput.h"

// This table should match order in EntapGlobals.h ENTAP_HEADERS enum
QueryData::EntapHeader QueryData::ENTAP_HEADER_INFO[] = {
        {"Unused",              false},                         // 0
        {"Query Sequence",      true},

        /* Frame Selection */
        {"Frame",               false},

        /* Expression Filtering */
        {"FPKM",                false},
        {"TPM",                 false},

        /* Similarity Search - General */
        {"Subject Sequence",    false},
        {"Percent Identical",   false},                          // 5
        {"Alignment Length",    false},
        {"Mismatches",          false},
        {"Gap Openings",        false},
        {"Query Start",         false},
        {"Query End",           false},                          // 10
        {"Subject Start",       false},
        {"Subject End",         false},
        {"E Value",             false},
        {"Coverage",            false},
        {"Description",         false},                          // 15
        {"Species",             false},
        {"Taxonomic Lineage",   false},
        {"Origin Database",     false},
        {"Contaminant",         false},
        {"Informative",         false},

        /* Similarity Search - UniProt */
        {"UniProt Database Cross Reference",        false},      // 20
        {"UniProt Additional Information",          false},
        {"UniProt KEGG Terms",                      false},
        {"UniProt GO Biological",                   false},
        {"UniProt GO Cellular",                     false},
        {"UniProt GO Molecular",                    false},

        /* Ontology - EggNOG */
        {"EggNOG Seed Ortholog",                   false},
        {"EggNOG Seed E-Value",                    false},
        {"EggNOG Seed Score",                      false},
        {"EggNOG Predicted Gene",                  false},
        {"EggNOG Tax Scope",                       false},
        {"EggNOG Tax Scope Max",                   false},
        {"EggNOG Member OGs",                      false},
        {"EggNOG Description",                     false},
        {"EggNOG BIGG Reaction",                   false},
        {"EggNOG KEGG Terms",                      false},
        {"EggNOG GO Biological",                   false},
        {"EggNOG GO Cellular",                     false},
        {"EggNOG GO Molecular" ,                   false},
        {"EggNOG Protein Domains",                 false},

        /* Ontology - InterProScan */
        {"IPScan GO Biological",                    false},     // 40
        {"IPScan GO Cellular",                      false},
        {"IPScan GO Molecular",                     false},
        {"IPScan Pathways",                         false},
        {"IPScan InterPro ID",                      false},
        {"IPScan Protein Database",                 false},
        {"IPScan Protein Description",              false},
        {"IPScan E-Value",                          false},

        /* Ontology - BUSCO */
        {"BUSCO ID",                                false},
        {"BUSCO Status",                            false},
        {"BUSCO Length",                            false},    // 50
        {"BUSCO Score",                             false},


        {"Unused",                                  false}
};


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
    uint16                                   shortest_len=0xFFFF;
    uint16                                   longest_len=0;
    uint16                                   len;
    fp64                                     avg_len;
    std::vector<uint16>                      sequence_lengths;
    std::pair<uint16, uint16>                n_vals;
    bool                                     is_complete;

    mTotalSequences = 0;
    mPipelineFlags  = 0;
    mDataFlags      = 0;
    mpSequences      = new QUERY_MAP_T;

    mpUserInput  = userinput;
    mpFileSystem = filesystem;

    mNoTrim        = mpUserInput->has_input(INPUT_FLAG_NO_TRIM);
    is_complete    = mpUserInput->has_input(INPUT_FLAG_COMPLETE);

    if (!mpFileSystem->file_exists(input_file)) {
        throw ExceptionHandler("Input transcriptome not found at: " + input_file,ERR_ENTAP_INPUT_PARSE);
    }

    out_name     = mpFileSystem->get_filename(input_file, true);
    out_new_path = PATHS(out_path,out_name);
    mpFileSystem->delete_file(out_new_path);

    // Set IS_PROTEIN flag if protein sequences found
    set_input_type(input_file);
    DATA_FLAG_GET(IS_PROTEIN) ? transcript_type = PROTEIN_FLAG : transcript_type = NUCLEO_FLAG;

    std::ifstream in_file(input_file);
    std::ofstream out_file(out_new_path,std::ios::out | std::ios::app);

    while (true) {
        std::getline(in_file, line);
        if (line.empty() && !in_file.eof()) continue;
        if (line.find(FileSystem::FASTA_FLAG) == 0 || in_file.eof()) {
            if (!seq_id.empty()) {
                if (in_file.eof()) {
                    out_file << line << std::endl;
                    sequence += line + "\n";
                }
                QuerySequence *query_seq = new QuerySequence(DATA_FLAG_GET(IS_PROTEIN),sequence, seq_id);
                if (is_complete) query_seq->setFrame(COMPLETE_FLAG);
                // Check for duplicate headers in input transcriptomes
                if (mpSequences->find(seq_id) != mpSequences->end()) {
                    throw ExceptionHandler("Duplicate headers in your input transcriptome: " + seq_id,
                        ERR_ENTAP_INPUT_PARSE);
                }
                mpSequences->emplace(seq_id, query_seq);
                count_seqs++;
                len = (uint16) query_seq->get_sequence_length();
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
    mTotalSequences = count_seqs;
    DATA_FLAG_GET(IS_PROTEIN)  ? mProteinLengthStart = total_len : mNucleoLengthStart = total_len;
    // first - n50, second - n90
    n_vals = calculate_N_vals(sequence_lengths, total_len);

    mpFileSystem->format_stat_stream(out_msg, "Transcriptome Statistics");
    out_msg <<
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
    mpFileSystem->print_stats(msg);
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
    FS_dprint("Transcriptome Lines - START");
    while(std::getline(in_file,line)) {
        if (line.empty()) continue;
        if (line_count++ > LINE_COUNT) break;
        if (line_count < SEQ_DPRINT_CONUT) FS_dprint(line);
        line.pop_back(); // Account for newline/other
        if (line.find('>') == std::string::npos) {
            for (char &c : line) {
                if (std::find(NUCLEO_MAP.begin(), NUCLEO_MAP.end(),toupper(c)) == NUCLEO_MAP.end())
                    deviations++;
            }
        }
    }
    FS_dprint("Transcriptome Lines - END");
    if (deviations > NUCLEO_DEV) {
        DATA_FLAG_SET(IS_PROTEIN);
    }
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

/**
 * ======================================================================
 * Function final_statistics(std::string &outpath)
 *
 * Description          - Calculates final statistical information after
 *                        completed execution
 *                      - Compiles stats on each stage of pipeline
 *
 * Notes                - None
 *
 * @param outpath       - Absolute path to statistics output directory
 *
 * @return              - None
 *
 * =====================================================================
 */
void QueryData::final_statistics(std::string &outpath) {
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
    ent_input_multi_int_t  ontology_flags;

    ontology_flags = mpUserInput->get_user_input<ent_input_multi_int_t>(INPUT_FLAG_ONTOLOGY);

    out_unannotated_nucl_path = PATHS(outpath, OUT_UNANNOTATED_NUCL);
    out_unannotated_prot_path = PATHS(outpath, OUT_UNANNOTATED_PROT);
    out_annotated_nucl_path   = PATHS(outpath, OUT_ANNOTATED_NUCL);
    out_annotated_prot_path   = PATHS(outpath, OUT_ANNOTATED_PROT);

    // Re-write these files
    mpFileSystem->delete_file(out_unannotated_nucl_path);
    mpFileSystem->delete_file(out_unannotated_prot_path);
    mpFileSystem->delete_file(out_annotated_nucl_path);
    mpFileSystem->delete_file(out_annotated_prot_path);

    std::ofstream file_unannotated_nucl(out_unannotated_nucl_path, std::ios::out | std::ios::app);
    std::ofstream file_unannotated_prot(out_unannotated_prot_path, std::ios::out | std::ios::app);
    std::ofstream file_annotated_nucl(out_annotated_nucl_path, std::ios::out | std::ios::app);
    std::ofstream file_annotated_prot(out_annotated_prot_path, std::ios::out | std::ios::app);

    for (auto &pair : *mpSequences) {
        count_total_sequences++;
        is_exp_kept = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_EXPRESSION_KEPT);
        is_prot = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN);
        is_hit = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_BLAST_HIT);
        is_ontology = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_FAMILY_ASSIGNED); // TODO Fix for interpro
        is_one_go = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_FAMILY_ONE_GO);
        is_one_kegg = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_FAMILY_ONE_KEGG);

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

    mpFileSystem->format_stat_stream(ss, "Final Annotation Statistics");
    ss <<
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
                case ONT_EGGNOG_DMND:
                    ss <<
                       "\nGene Families"        <<
                       "\n\tTotal unique sequences with family assignment: "    << count_ontology   <<
                       "\n\tTotal unique sequences without family assignment: " << count_no_ontology<<
                       "\n\tTotal unique sequences with at least one GO term: " << count_one_go     <<
                       "\n\tTotal unique sequences with at least one pathway (KEGG) assignment: "   << count_one_kegg;
                    break;
                case ONT_INTERPRO_SCAN:
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
    mpFileSystem->print_stats(out_msg);
}

std::string QueryData::trim_sequence_header(std::string &header, std::string line) {
    std::string   sequence;
    int16         pos;

    if (line.find('>') != std::string::npos) {
        pos = (int16) line.find('>');
    } else pos = -1;
    if (mNoTrim) {
        // No trimming, just remove spaces
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        header = line.substr(pos+1);
        sequence = line + "\n";
    } else {
        // Trim to first space
        if (line.find(' ') != std::string::npos) {
            header = line.substr(pos+1, line.find(' ')-(pos+1));
        } else header = line.substr(pos+1);
        sequence = ">" + header + "\n";
    }
    return sequence;
}

QUERY_MAP_T* QueryData::get_sequences_ptr() {
    return this->mpSequences;
}

QueryData::~QueryData() {
    std::string path;
    FS_dprint("Killing Object - QueryData");
    for (auto &mpSequence : *mpSequences) {
        delete mpSequence.second;
        mpSequence.second = nullptr;
    }
    FS_dprint("QuerySequence data freed");
    delete mpSequences;

    // Cleanup files in case it was interrupted
    for (auto &pair : mAlignmentFiles) {
        path = pair.first;
        end_alignment_files(path);
    }
}

bool QueryData::DATA_FLAG_GET(DATA_FLAGS flag) {
    return (mDataFlags & flag) != 0;
}

void QueryData::DATA_FLAG_SET(DATA_FLAGS flag) {
    mDataFlags |= flag;
}

void QueryData::DATA_FLAG_CLEAR(DATA_FLAGS flag) {
    mDataFlags &= ~flag;
}

QuerySequence *QueryData::get_sequence(std::string &query_id) {
    QUERY_MAP_T::iterator it = mpSequences->find(query_id);
    if (it != mpSequences->end()) {
        // Sequence found, return
        return it->second;
    } else {
        // Sequence NOT found retun null
        return nullptr;
    }
}

bool QueryData::start_alignment_files(std::string &base_path, std::vector<ENTAP_HEADERS> &headers, uint8 lvl,
                                        std::vector<FileSystem::ENT_FILE_TYPES> &types) {
    bool ret;
    OutputFileData outputFileData = OutputFileData();
    outputFileData.headers = headers;
    outputFileData.go_level = lvl;
    outputFileData.file_types = types;

    base_path = base_path + "_lvl" + std::to_string(lvl);

    // add this path to map if it does not exist, otherwise skip
    if (mAlignmentFiles.find(base_path) == mAlignmentFiles.end()) {

        mAlignmentFiles.emplace(base_path, outputFileData);
        // Generate files for each data type
        for (FileSystem::ENT_FILE_TYPES type : types) {
            if ((type == FileSystem::ENT_FILE_FASTA_FAA || type == FileSystem::ENT_FILE_FASTA_FNN) && lvl != 0) {
                // Do NOT create files for anything other than 0 for FAA or FNN
                continue;
            } else {
                mAlignmentFiles.at(base_path).file_streams[type] =
                        new std::ofstream(base_path + mpFileSystem->get_extension(type), std::ios::out | std::ios::app);
                // Initialize headers or any other generic stuff
                initialize_file(mAlignmentFiles.at(base_path).file_streams[type], headers, type);
            }
        }
        ret = true;
    } else {
        // File base path already exists!
        ret = false;
    }
    return ret;
}

bool QueryData::end_alignment_files(std::string &base_path) {
    // Cleanup/close files

    for (std::ofstream* file_ptr : mAlignmentFiles.at(base_path).file_streams) {
        // some are unused such as 0
        if (file_ptr != nullptr) {
            file_ptr->close();
            delete file_ptr;
        }
    }
    mAlignmentFiles.erase(base_path);
    return true;
}

bool QueryData::add_alignment_data(std::string &base_path, QuerySequence *querySequence, QueryAlignment *alignment) {
    bool ret = false;

    // Cycle through output file types for this path
    for (FileSystem::ENT_FILE_TYPES type : mAlignmentFiles.at(base_path).file_types) {

        if (mAlignmentFiles.at(base_path).file_streams[type] == nullptr) continue;

        switch (type) {

            case FileSystem::ENT_FILE_DELIM_TSV:
                if (alignment == nullptr) {
                    *mAlignmentFiles.at(base_path).file_streams[type] <<
                        get_delim_data_sequence(mAlignmentFiles.at(base_path).headers,FileSystem::DELIM_TSV,
                        mAlignmentFiles.at(base_path).go_level, querySequence) << std::endl;
                } else {
                    *mAlignmentFiles.at(base_path).file_streams[type] <<
                        get_delim_data_alignment(mAlignmentFiles.at(base_path).headers,FileSystem::DELIM_TSV,
                        mAlignmentFiles.at(base_path).go_level, alignment) << std::endl;
                }
                break;

            case FileSystem::ENT_FILE_DELIM_CSV:
                if (alignment == nullptr) {
                    *mAlignmentFiles.at(base_path).file_streams[type] <<
                        get_delim_data_sequence(mAlignmentFiles.at(base_path).headers, FileSystem::DELIM_CSV,
                        mAlignmentFiles.at(base_path).go_level, querySequence) << std::endl;
                } else {
                    *mAlignmentFiles.at(base_path).file_streams[type] <<
                        get_delim_data_alignment(mAlignmentFiles.at(base_path).headers,FileSystem::DELIM_CSV,
                        mAlignmentFiles.at(base_path).go_level, alignment) << std::endl;
                }
                break;

            case FileSystem::ENT_FILE_FASTA_FAA:
                if (!querySequence->get_sequence_p().empty())
                    *mAlignmentFiles.at(base_path).file_streams[type] << querySequence->get_sequence_p() << std::endl;
                break;

            case FileSystem::ENT_FILE_FASTA_FNN:
                if (!querySequence->get_sequence_n().empty())
                    *mAlignmentFiles.at(base_path).file_streams[type] << querySequence->get_sequence_n() << std::endl;
                break;

            default:
                FS_dprint("ERROR unhandled file type (add_alignment_data): " + std::to_string(type));
                break;
        }
    }
    return ret;
}

bool QueryData::is_protein_data() {
    return DATA_FLAG_GET(IS_PROTEIN);
}

void QueryData::set_is_protein_data(bool val) {
    DATA_FLAG_CHANGE(IS_PROTEIN, val);
}

void QueryData::set_is_success_frame_selection(bool val) {
    DATA_FLAG_CHANGE(SUCCESS_FRAME_SEL, val);
}

void QueryData::set_is_success_expression(bool val) {
    DATA_FLAG_CHANGE(SUCCESS_EXPRESSION, val);
}

void QueryData::set_is_success_sim_search(bool val) {
    DATA_FLAG_CHANGE(SUCCESS_SIM_SEARCH, val);
}

void QueryData::set_is_success_ontology(bool val) {
    DATA_FLAG_CHANGE(SUCCESS_ONTOLOGY, val);
}

void QueryData::DATA_FLAG_CHANGE(QueryData::DATA_FLAGS flag, bool val) {
    if (val) {
        DATA_FLAG_SET(flag);
    } else {
        DATA_FLAG_CLEAR(flag);
    }
}

void QueryData::set_is_uniprot(bool val) {
    DATA_FLAG_CHANGE(UNIPROT_MATCH, val);
}

void QueryData::header_set(ENTAP_HEADERS header, bool val) {
    ENTAP_HEADER_INFO[header].print_header = val;
}

/**
 * ======================================================================
 * Function std::string QueryData::get_delim_data_sequence(std::vector<ENTAP_HEADERS> &headers,
 *                                              char delim, uint8 lvl,
 *                                              QuerySequence *sequence)
 *
 * Description          - Accesses QuerySequence object to pull all relevant
 *                        Header data for input headers
 *
 * Notes                - None
 *
 * @param headers       - EnTAP headers we would like to pull data for
 * @param delim         - Character to use as deliminator for return string
 * @param lvl           - Gene Ontology level we would like data for
 * @param sequence      - Pointer to QuerySequence object to pull data from
 *
 * @return              - String of relevant header data with specified deliminator
 *
 * =====================================================================
 */
std::string QueryData::get_delim_data_sequence(std::vector<ENTAP_HEADERS> &headers, char delim,
                                               uint8 lvl, QuerySequence *sequence) {
    std::string val;
    std::stringstream ret;

    for (ENTAP_HEADERS header : headers) {
        if (ENTAP_HEADER_INFO[header].print_header) {
            sequence->get_header_data(val, header, lvl);
            ret << val << delim;
        }
    }
    return ret.str();
}

/**
 * ======================================================================
 * Function std::string QueryData::get_delim_data_alignment(std::vector<ENTAP_HEADERS> &headers,
 *                                              char delim, uint8 lvl,
 *                                              QueryAlignment *alignment)
 *
 * Description          - Accesses QueryAlignment object to pull all relevant
 *                        Header data for input headers
 *
 * Notes                - None
 *
 * @param headers       - EnTAP headers we would like to pull data for
 * @param delim         - Character to use as deliminator for return string
 * @param lvl           - Gene Ontology level we would like data for
 * @param alignment     - Pointer to QueryAlignment object to pull data from
 *
 * @return              - String of relevant header data with specified delminator
 *
 * =====================================================================
 */
std::string QueryData::get_delim_data_alignment(std::vector<ENTAP_HEADERS> &headers, char delim,
                                                uint8 lvl, QueryAlignment *alignment) {
    std::string val;
    std::stringstream ret;

    for (ENTAP_HEADERS header : headers) {
        if (ENTAP_HEADER_INFO[header].print_header) {
            alignment->get_header_data(header, val, lvl);
            ret << val << delim;
        }
    }
    return ret.str();
}

/**
 * ======================================================================
 * Function bool QueryData::initialize_file(std::ofstream *file_stream,
 *                               std::vector<ENTAP_HEADERS> &headers,
 *                               FileSystem::ENT_FILE_TYPES type)
 *
 * Description          - Initializes output files depending on their type
 *                        (ex: deliminated files will print headers as the
 *                        first line of the file)
 *
 * Notes                - None
 *
 * @param file_stream   - File stream to print to
 * @param headers       - EnTAP headers we are considering for this file (delim)
 * @param type          - Type of file we want to initialize
 *
 * @return              - None
 *
 * =====================================================================
 */
bool QueryData::initialize_file(std::ofstream *file_stream, std::vector<ENTAP_HEADERS> &headers,
                                 FileSystem::ENT_FILE_TYPES type) {
    bool ret=true;

    switch (type) {
        case FileSystem::ENT_FILE_DELIM_TSV:
            for (ENTAP_HEADERS &header: headers) {
                if (ENTAP_HEADER_INFO[header].print_header) {
                    *file_stream << ENTAP_HEADER_INFO[header].title << FileSystem::DELIM_TSV;
                }
            }
            *file_stream << std::endl;
            break;

        case FileSystem::ENT_FILE_DELIM_CSV:
            for (ENTAP_HEADERS &header: headers) {
                if (ENTAP_HEADER_INFO[header].print_header) {
                    *file_stream << ENTAP_HEADER_INFO[header].title << FileSystem::DELIM_CSV;
                }
            }
            *file_stream << std::endl;
            break;

            // Fasta files do not need initialization
        case FileSystem::ENT_FILE_FASTA_FAA:
        case FileSystem::ENT_FILE_FASTA_FNN:
            break;

        default:
            FS_dprint("ERROR unhandled file type (initialize file): " + std::to_string(type));
            ret = false;
            break;
    }
    return ret;
}

QUERY_MAP_T QueryData::get_specific_sequences(uint32 flags) {
    QUERY_MAP_T ret_map = QUERY_MAP_T();

    for (auto &pair : *mpSequences) {
        // Check if this sequence has any of the flags we want
        if (pair.second->QUERY_FLAG_CONTAINS(flags)){
            // Yes, add to map
            ret_map.emplace(pair.first, pair.second);
        }
    }
    return ret_map;
}
