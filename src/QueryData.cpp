/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2023, Alexander Hart, Dr. Jill Wegrzyn
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
        {"FPKM",                        false},
        {"TPM",                         false},
        {"Expression Effective Length", false},

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
 * Function QueryData::QueryData(std::string &input_path,
 *                               UserInput *userInput,
 *                               FileSystem *fileSystem)
 *
 * Description          - Parses input transcriptome and converts to map of
 *                        each query sequence
 *                      - This map is passed throughout EnTAP execution and
 *                        updated throughout
 *
 * Notes                - Used during Config, will NOT print output
 *
 * @param input_path    - Absolute path to input transcriptome
 * @param userInput     - Pointer to UserInput object
 * @param fileSystem    - Pointer to FileSystem object
 *
 * @return              - None
 *
 * =====================================================================
 */
QueryData::QueryData(std::string &input_path, UserInput *userInput, FileSystem *fileSystem) {
    FS_dprint("Spawn Object - QueryData");

    init_params(fileSystem, userInput); // Initialize member variables
    try {
        generate_transcriptome(input_path, false, "");
    } catch (const ExceptionHandler &e) {
        throw e;
    }
}


/**
 * ======================================================================
 * Function QueryData::QueryData(std::string &input_file,
 *                              std::string &out_path,
 *                              userInput  *userinput,
 *                              FileSystem *filesystem,
 *
 * Description          - Parses input transcriptome and converts to map of
 *                        each query sequence
 *                      - This map is passed throughout EnTAP execution and
 *                        updated throughout
 *
 * Notes                - Will print transcriptome to output file + generate log
 *
 * @param input_file    - Path to input transcriptome
 * @param out_path      - Path to output formatted transcriptome
 * @param userinput     - Pointer to UserInput object
 * @param filesystem    - Pointer to FileSystem object
 *
 * @return              - Sets input_file to out_path
 *
 * =====================================================================
 */
QueryData::QueryData(std::string &input_file, std::string &out_path, UserInput *userinput,
                    FileSystem* filesystem) {
    FS_dprint("Spawn Object - QueryData");

    init_params(filesystem, userinput); // Initialize member variables

    try{
        generate_transcriptome(input_file, true, out_path);
    } catch (const ExceptionHandler &e) {
        throw e;
    }

    input_file = out_path;
}


void QueryData::generate_transcriptome(std::string &input_path, bool print_output, std::string output_path) {

    std::string                     line;                   // Line read from input trancriptome (temp)
    std::string                     seq_id;                 // Sequence ID of each sequence in input (temp)
    std::string                     sequence;               // Sequence of each in input transcriptome (temp)
    std::string                     longest_seq;            // Sequence ID of longest sequence
    std::string                     shortest_seq;           // Sequence ID of shortest sequence
    std::ifstream                   in_file;                // INPUT stream of input file
    std::ofstream                   out_file;               // OUTPUT stream of output file
    std::stringstream               out_msg;                // Output string stream to log file
    uint32                          count_seqs=0;           // Total number of sequences
    uint64                          total_len=0;            // Total length (BP) of transcriptome
    uint16                          shortest_len=0xFFFF;    // Length (BP) of shortest sequence
    uint16                          longest_len=0;          // Length (BP) of longest sequence
    uint16                          len;
    fp64                            avg_len;
    std::vector<uint16>             sequence_lengths;
    std::pair<uint16, uint16>       n_vals;

    FS_dprint("Generating transcriptome mappings...");

    // Verify input
    if (!mpFileSystem->file_exists(input_path)) {
        throw ExceptionHandler("Input transcriptome not found at: " + input_path,ERR_ENTAP_INPUT_PARSE);
    } else if (mpFileSystem->file_empty(input_path)) {
        throw ExceptionHandler("Input transcriptome is empty at: " + input_path, ERR_ENTAP_INPUT_PARSE);
    }

    // Set IS_PROTEIN flag if protein sequences detected
    set_input_type(input_path);

    // Ensure we can open input file
    try {
        in_file.open(input_path);  // Open transcriptome
    } catch(std::exception &e) {
        FS_dprint("ERROR opening input transcriptome.: " + std::string(e.what()));
        throw ExceptionHandler("Error opening input transcriptome, check permissions at: " +
            input_path, ERR_ENTAP_INPUT_PARSE);
    }

    // Ensure we can write to a file
    if (print_output) {
        try {
            mpFileSystem->delete_file(output_path);
            out_file.open(output_path,std::ios::out | std::ios::app);
        } catch (std::exception &e) {
            throw ExceptionHandler("ERROR writing to output file, check permissions at: " + output_path,
                                   ERR_ENTAP_INPUT_PARSE);
        }
    }

    // Being to parse transcriptome and generate data map
    while (true) {
        std::getline(in_file, line);
        if (line.empty() && !in_file.eof()) continue;
        if (line.find(FileSystem::FASTA_FLAG) == 0 || in_file.eof()) {
            if (!seq_id.empty()) {
                if (in_file.eof()) {
                    if (print_output) out_file << line << std::endl;
                    sequence += line + "\n";
                }
                QuerySequence *query_seq = new QuerySequence(DATA_FLAG_GET(IS_PROTEIN),sequence, seq_id);
                if (mIsComplete) query_seq->setFrame(COMPLETE_FLAG);
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
                    longest_len = len;
                    longest_seq = seq_id;
                }
                if (len < shortest_len) {
                    shortest_len = len;
                    shortest_seq = seq_id;
                }
                sequence_lengths.push_back(len);
            }
            if (in_file.eof()) break;
            sequence = trim_sequence_header(seq_id, line);
            if (print_output) out_file << sequence;
        } else {
            if (print_output) out_file << line << std::endl;
            sequence += line + "\n";
        }
    }
    // Close files
    in_file.close();
    if (print_output) out_file.close();

    avg_len = total_len / count_seqs;
    mTotalSequences = count_seqs;
    DATA_FLAG_GET(IS_PROTEIN)  ? mProteinLengthStart = total_len : mNucleoLengthStart = total_len;
    // first - n50, second - n90
    n_vals = calculate_N_vals(sequence_lengths, total_len);

    // Print output to Log file and debug
    if (print_output) {
        mpFileSystem->format_stat_stream(out_msg, "Transcriptome Statistics");
        out_msg <<
                mTranscriptTypeStr << " sequences found"          <<
                "\nTotal sequences: "                          << count_seqs    <<
                "\nTotal length of transcriptome(bp): "        << total_len     <<
                "\nAverage sequence length(bp): "              << avg_len       <<
                "\nn50: "                                      << n_vals.first  <<
                "\nn90: "                                      << n_vals.second <<
                "\nLongest sequence(bp): " << longest_len << " ("<<longest_seq<<")"<<
                "\nShortest sequence(bp): "<< shortest_len<<" ("<<shortest_seq<<")";
        if (mIsComplete)out_msg<<"\nAll sequences ("<<count_seqs<<") were flagged as complete genes";
        std::string msg = out_msg.str();
        mpFileSystem->print_stats(msg);
    }
    FS_dprint("Success!");
}


void QueryData::init_params(FileSystem *fileSystem, UserInput *userInput) {
    mTotalSequences = 0;
    mDataFlags      = 0;
    mpSequences      = new QUERY_MAP_T;

    mpUserInput  = userInput;
    mpFileSystem = fileSystem;

    mNoTrim        = mpUserInput->has_input(INPUT_FLAG_NO_TRIM);
    mIsComplete    = mpUserInput->has_input(INPUT_FLAG_COMPLETE);
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
        mTranscriptTypeStr = PROTEIN_FLAG;
    } else{
        DATA_FLAG_SET(IS_NUCLEOTIDE);
        mTranscriptTypeStr = NUCLEO_FLAG;
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
    uint32                 count_total_kept_sequences=0;
    uint32                 count_exp_kept=0;
    uint32                 count_exp_reject=0;
    uint32                 count_frame_kept=0;
    uint32                 count_frame_rejected=0;
    uint32                 count_sim_hits=0;
    uint32                 count_sim_no_hits=0;
    uint32                 count_sim_contam=0;          // Unique alignments flagged as a contaminant
    uint32                 count_ontology=0;
    uint32                 count_no_ontology=0;
    uint32                 count_one_go=0;
    uint32                 count_one_kegg=0;
    uint32                 count_sim_only=0;
    uint32                 count_ontology_only=0;
    uint32                 count_TOTAL_ann=0;
    uint32                 count_TOTAL_unann=0;
    uint32                 count_TOTAL_unann_kept=0;    // Total sequences unannotated if kept after expression/frame selection
    std::string            out_unannotated_nucl_path;
    std::string            out_unannotated_prot_path;
    std::string            out_annotated_nucl_path;
    std::string            out_annotated_prot_path;
    std::string            out_msg;
    bool                   is_exp_kept;
    bool                   is_prot;
    bool                   is_frame_kept;
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
        is_frame_kept = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_FRAME_KEPT);
        is_prot = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN);
        is_hit = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_BLAST_HIT);
        is_ontology = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_FAMILY_ASSIGNED); // TODO Fix for interpro
        is_one_go = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_FAMILY_ONE_GO);
        is_one_kegg = pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_FAMILY_ONE_KEGG);

        if (pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_CONTAMINANT)) {
            count_sim_contam++; // This is assumed as retained and currently only from sim search
        }

        is_exp_kept ? count_exp_kept++ : count_exp_reject++;
        is_prot ? count_frame_kept++ : count_frame_rejected++;
        is_hit ? count_sim_hits++ : count_sim_no_hits++;
        is_ontology ? count_ontology++ : count_no_ontology++;
        if (is_one_go) count_one_go++;
        if (is_one_kegg) count_one_kegg++;

        if (is_hit && !is_ontology) count_sim_only++;
        if (!is_hit && is_ontology) count_ontology_only++;

        if (is_exp_kept && is_frame_kept) {
            count_total_kept_sequences++;
        }

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
            if (is_exp_kept && is_frame_kept) {
                count_TOTAL_unann_kept++;
            }
            else {
                count_TOTAL_unann++;
            }
        }
    }

    file_unannotated_nucl.close();
    file_unannotated_prot.close();
    file_annotated_nucl.close();
    file_annotated_prot.close();

    mpFileSystem->format_stat_stream(ss, "Final Annotation Statistics");
    ss <<
       "Total Input Sequences: "                  << count_total_sequences;

    if (DATA_FLAG_GET(SUCCESS_EXPRESSION)) {
        ss <<
           "\nExpression Analysis" <<
           "\n\tRetained sequences: "  << count_exp_kept    << " (" <<
                (((fp64) count_exp_kept / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)" <<
           "\n\tLost sequences: "  << count_exp_reject << " (" <<
                (((fp64) count_exp_reject / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)";
    }
    if (DATA_FLAG_GET(SUCCESS_FRAME_SEL)) {
        ss <<
           "\nFrame Selection"              <<
           "\n\tTotal sequences retained: " << count_frame_kept     << " (" <<
                (((fp64) count_frame_kept / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)" <<
           "\n\tTotal sequences removed: "  << count_frame_rejected << " (" <<
                (((fp64) count_frame_rejected / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)";
    }
    if (DATA_FLAG_GET(SUCCESS_SIM_SEARCH)) {
        ss <<
           "\nSimilarity Search"                               <<
           "\n\tTotal unique sequences with an alignment: "    << count_sim_hits   << " (" <<
                 (((fp64) count_sim_hits / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)" <<
           "\n\t\tTotal alignments flagged as a contaminant: "   << count_sim_contam << " (" <<
                 (((fp64) count_sim_contam / count_sim_hits) * ENTAP_PERCENT) << "% of total unique alignments)" <<
           "\n\t\tTotal alignments NOT flagged as a contaminant: "   << (count_sim_hits - count_sim_contam) << " (" <<
                 (((fp64) (count_sim_hits - count_sim_contam) / count_sim_hits) * ENTAP_PERCENT) << "% of total unique alignments)" <<
           "\n\tTotal unique sequences without an alignment: " << count_sim_no_hits << " (" <<
                 (((fp64) count_sim_no_hits / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)";
    }
    if (DATA_FLAG_GET(SUCCESS_ONTOLOGY)) {
        for (uint16 flag : ontology_flags) {
            switch (flag) {
                case ONT_EGGNOG_DMND:
                    ss <<
                       "\nGene Families"        <<
                       "\n\tTotal unique sequences with family assignment: "    << count_ontology   << " (" <<
                            (((fp64) count_ontology / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)" <<
                       "\n\tTotal unique sequences without family assignment: " << count_no_ontology<< " (" <<
                            (((fp64) count_no_ontology / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)" <<
                       "\n\tTotal unique sequences with at least one GO term: " << count_one_go     << " (" <<
                            (((fp64) count_one_go / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)" <<
                       "\n\tTotal unique sequences with at least one pathway (KEGG) assignment: "   << count_one_kegg << " (" <<
                           (((fp64) count_one_kegg / count_total_sequences) * ENTAP_PERCENT) << "% of total input sequences)";
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
       "\n\tTotal Retained Sequences (After Filtering And/Or Frame Selection): " << count_total_kept_sequences <<
       "\n\tTotal unique sequences annotated (similarity search alignments only): "      << count_sim_only      << " (" <<
            (((fp64) count_sim_only / count_total_kept_sequences) * ENTAP_PERCENT) << "% of total retained)" <<
       "\n\tTotal unique sequences annotated (gene family assignment only): "            << count_ontology_only << " (" <<
            (((fp64) count_ontology_only / count_total_kept_sequences) * ENTAP_PERCENT) << "% of total retained)" <<
       "\n\tTotal unique sequences annotated (gene family and/or similarity search): "   << count_TOTAL_ann     << " (" <<
            (((fp64) count_TOTAL_ann / count_total_kept_sequences) * ENTAP_PERCENT) << "% of total retained)" <<
       "\n\t\tTotal alignments flagged as a contaminant (gene family and/or similarity search): " << count_sim_contam << " (" <<
            (((fp64) count_sim_contam / count_sim_hits) * ENTAP_PERCENT) << "% of total unique alignments)" <<
       "\n\t\tTotal alignments NOT flagged as a contaminant (gene family and/or similarity search): " << (count_sim_hits - count_sim_contam) << " (" <<
            (((fp64) (count_sim_hits - count_sim_contam) / count_sim_hits) * ENTAP_PERCENT) << "% of total unique alignments)" <<
       "\n\tTotal unique sequences unannotated (gene family and/or similarity search): " << count_TOTAL_unann_kept << " ("
            << (((fp64) count_TOTAL_unann_kept / count_total_kept_sequences) * ENTAP_PERCENT) << "% of total retained)";

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
        std::string::iterator it = std::find_if(line.begin(), line.end(), isspace);
        if (it != line.end()) {
            header = line.substr(pos+1, std::distance(line.begin(), it)-(pos+1));
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
    if (mpSequences != nullptr) {
        for (auto &mpSequence : *mpSequences) {
            delete mpSequence.second;
            mpSequence.second = nullptr;
        }
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

bool QueryData::start_alignment_files(const std::string &base_path, const std::vector<ENTAP_HEADERS> &headers, const vect_uint16_t &go_levels,
                                      const std::vector<FileSystem::ENT_FILE_TYPES> &types) {
    bool ret;
    bool is_fasta_generated;
    std::string final_file_path; // Final path may change depending on go levels, etc...

    OutputFileData outputFileData = OutputFileData();
    outputFileData.headers = headers;
    outputFileData.go_levels = go_levels;
    outputFileData.file_types = types;

    // add this path to map if it does not exist, otherwise skip
    if (mAlignmentFiles.find(base_path) == mAlignmentFiles.end()) {
        mAlignmentFiles.emplace(base_path, outputFileData);
        // Generate files for each data type
        for (FileSystem::ENT_FILE_TYPES type : types) {
            // Generate files for each go level
            is_fasta_generated = false; // FASTA only has one file across go levels

            // Just skip out of range for now (checked when verifying user input)
            if (type >= FileSystem::ENT_FILE_OUTPUT_FORMAT_MAX) continue;

            for (uint16 level : go_levels) {

                // Just skip out of range for now (checked when verifying user input)
                if (level >= UserInput::MAX_GO_LEVEL) continue;

                switch (type) {

                    // For TSV,CSV we want to create a new file for every go level
                    case FileSystem::ENT_FILE_DELIM_TSV:
                    case FileSystem::ENT_FILE_DELIM_CSV:
                        final_file_path = base_path + APPEND_GO_LEVEL_STR +
                                std::to_string(level) + mpFileSystem->get_extension(type);
                        mAlignmentFiles.at(base_path).file_streams[type][level] =
                                new std::ofstream(final_file_path, std::ios::out | std::ios::app);
                        initialize_file(mAlignmentFiles.at(base_path).file_streams[type][level],
                                        outputFileData.headers, type);
                        break;

                    // For FASTA, we do not want to create multiple for each go level
                    case FileSystem::ENT_FILE_FASTA_FAA:
                    case FileSystem::ENT_FILE_FASTA_FNN:
                        if (!is_fasta_generated) {
                            final_file_path = base_path + mpFileSystem->get_extension(type);
                            is_fasta_generated = true;
                            mAlignmentFiles.at(base_path).file_streams[type][DEFAULT_GO_LEVEL] =
                                    new std::ofstream(final_file_path, std::ios::out | std::ios::app);
                            initialize_file(mAlignmentFiles.at(base_path).file_streams[type][DEFAULT_GO_LEVEL],
                                            outputFileData.headers, type);
                        }
                        break;

                    case FileSystem::ENT_FILE_GENE_ENRICH_EFF_LEN:
                        final_file_path = base_path + APPEND_GO_LEVEL_STR +
                                          std::to_string(level) + APPEND_ENRICH_GENE_ID_LEN +
                                          mpFileSystem->get_extension(type);
                        mAlignmentFiles.at(base_path).file_streams[type][level] =
                                new std::ofstream(final_file_path, std::ios::out | std::ios::app);
                        initialize_file(mAlignmentFiles.at(base_path).file_streams[type][level],
                                        outputFileData.headers, type);
                        break;

                    case FileSystem::ENT_FILE_GENE_ENRICH_GO_TERM:
                        final_file_path = base_path + APPEND_GO_LEVEL_STR +
                                          std::to_string(level) + APPEND_ENRICH_GENE_ID_GO +
                                          mpFileSystem->get_extension(type);
                        mAlignmentFiles.at(base_path).file_streams[type][level] =
                                new std::ofstream(final_file_path, std::ios::out | std::ios::app);
                        initialize_file(mAlignmentFiles.at(base_path).file_streams[type][level],
                                        outputFileData.headers, type);
                        break;

                    default:
                        break;
                }
            }
        }
        ret = true;
        mAlignmentFiles.emplace(base_path, outputFileData);
        FS_dprint("Alignment file started: " + base_path);
    } else {
        // File base path already exists!
        FS_dprint("Alignment file skipped, already exists: " + base_path);
        ret = false;
    }
    return ret;
}

bool QueryData::end_alignment_files(std::string &base_path) {
    // Cleanup/close files

    auto file_data = mAlignmentFiles.find(base_path);

    if (file_data == mAlignmentFiles.end()) return true;

    for (FileSystem::ENT_FILE_TYPES type : file_data->second.file_types) {
        for (std::ofstream* file_ptr : mAlignmentFiles.at(base_path).file_streams[type]) {
            // some are unused such as 0
            if (file_ptr != nullptr) {
                file_ptr->close();
                SAFE_DELETE(file_ptr);
            }
        }
    }

    mAlignmentFiles.erase(base_path);
    FS_dprint("Alignment file ended: " + base_path);
    return true;
}

bool QueryData::add_alignment_data(std::string &base_path, QuerySequence *querySequence, QueryAlignment *alignment) {
    bool ret = false;
    bool ignore_go_levels = false; // Only want the query sequence and ignore multiple go levels to prevent duplication

    auto file_data = mAlignmentFiles.find(base_path);

    if (file_data == mAlignmentFiles.end()) return ret; // EXIT ERROR alignment files missing

    // Cycle through output file types for this path
    for (FileSystem::ENT_FILE_TYPES type : file_data->second.file_types) {

        if (file_data->second.file_streams[type] == nullptr) continue;

        // loop through go levels
        for (uint16 go_level : file_data->second.go_levels) {
            switch (type) {

                case FileSystem::ENT_FILE_DELIM_TSV:
                    if (file_data->second.file_streams[type][go_level] == nullptr) continue;
                    if (alignment == nullptr) {
                        *mAlignmentFiles.at(base_path).file_streams[type][go_level] <<
                                                get_delim_data_sequence(file_data->second.headers,FileSystem::DELIM_TSV,
                                                go_level, querySequence) << std::endl;
                    } else {
                        *mAlignmentFiles.at(base_path).file_streams[type][go_level] <<
                                                get_delim_data_alignment(mAlignmentFiles.at(base_path).headers,FileSystem::DELIM_TSV,
                                                go_level, alignment) << std::endl;
                    }
                    break;

                case FileSystem::ENT_FILE_DELIM_CSV:
                    if (file_data->second.file_streams[type][go_level] == nullptr) continue;
                    if (alignment == nullptr) {
                        *mAlignmentFiles.at(base_path).file_streams[type][go_level] <<
                                                get_delim_data_sequence(file_data->second.headers, FileSystem::DELIM_CSV,
                                                go_level, querySequence) << std::endl;
                    } else {
                        *mAlignmentFiles.at(base_path).file_streams[type][go_level] <<
                                                get_delim_data_alignment(file_data->second.headers,FileSystem::DELIM_CSV,
                                                go_level, alignment) << std::endl;
                    }
                    break;

                // Skip FASTA files if QuerySequence has already been added for alignment?
                case FileSystem::ENT_FILE_FASTA_FAA:
                    if (file_data->second.file_streams[type][DEFAULT_GO_LEVEL] == nullptr) continue;
                    if (!querySequence->get_sequence_p().empty())
                        *mAlignmentFiles.at(base_path).file_streams[type][DEFAULT_GO_LEVEL] << querySequence->get_sequence_p() << std::endl;
                    ignore_go_levels = true;
                    break;

                case FileSystem::ENT_FILE_FASTA_FNN:
                    if (file_data->second.file_streams[type][DEFAULT_GO_LEVEL] == nullptr) continue;
                    if (!querySequence->get_sequence_n().empty())
                        *mAlignmentFiles.at(base_path).file_streams[type][DEFAULT_GO_LEVEL] << querySequence->get_sequence_n() << std::endl;
                    ignore_go_levels = true;
                    break;

                /*
                 * For Gene Enrichment with GO Terms, we want gene ID and Gene Ontology term
                 *  on the same line, tab deliminated
                 *
                 * If there are multiple GO terms to a gene id, there will be a new line for each GO
                 *  term.
                 *
                 * Example:
                 *  Gene_01 GO:0921901
                 *  Gene_01 GO:4321432
                 *  Gene_01 GO:4433233
                 *  Gene_02 GO:0000001
                 *
                 * */
                case FileSystem::ENT_FILE_GENE_ENRICH_GO_TERM:
                {
                    if (file_data->second.file_streams[type][go_level] == nullptr) continue;
                    if (!querySequence->is_kept()) continue;

                    go_format_t go_terms = querySequence->get_go_terms();
                    if (!go_terms.empty()) {
                        for (GoEntry const &entry : go_terms) {
                            if (!entry.go_id.empty()) {
                                if ((entry.level_int >= go_level && entry.level_int != GoEntry::UNKNOWN_LVL) ||
                                    go_level == 0) {
                                    *mAlignmentFiles.at(base_path).file_streams[type][go_level] <<
                                            querySequence->getMSequenceID() <<
                                            FileSystem::DELIM_TSV <<
                                            entry.go_id <<
                                            std::endl;
                                }
                            }
                        }
                    }
                    break;
                }

                case FileSystem::ENT_FILE_GENE_ENRICH_EFF_LEN:
                {
                    if (file_data->second.file_streams[type][go_level] == nullptr) continue;
                    if (!querySequence->is_kept()) continue;

                    if (querySequence->contains_go_level(go_level) || querySequence->contains_go_level(GoEntry::UNKNOWN_LVL)) {
                        *mAlignmentFiles.at(base_path).file_streams[type][go_level] <<
                                querySequence->getMSequenceID() <<
                                FileSystem::DELIM_TSV <<
                                querySequence->getMEffectiveLength() <<
                                std::endl;
                    }
                    break;
                }

                default:
                    FS_dprint("ERROR unhandled file type (add_alignment_data): " + std::to_string(type));
                    break;
            }

            // Certain file types do not care about multiple GO levels so we are going to ignore duplicate printing
            //  will restructure file printing later to not have this issue
            if (ignore_go_levels) {
                break;
            }
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
            // TO support tidyverse format we want empty data in TSV format to print 'NA'
            if (val.empty() && delim == FileSystem::DELIM_TSV) {
                    ret << FileSystem::TIDYVERSE_TSV_NULL << delim;
            }
            else {
                ret << val << delim;
            }
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
 * @return              - False if error
 *
 * =====================================================================
 */
bool QueryData::initialize_file(std::ofstream *file_stream, std::vector<ENTAP_HEADERS> &headers,
                                 FileSystem::ENT_FILE_TYPES type) {
    bool ret=true;
    if (file_stream == nullptr) return false;

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

        case FileSystem::ENT_FILE_GENE_ENRICH_EFF_LEN:
            *file_stream << HEADER_ENRICH_GENE_ID << FileSystem::DELIM_TSV
                         << HEADER_ENRICH_LENGTH << std::endl;
            break;

        case FileSystem::ENT_FILE_GENE_ENRICH_GO_TERM:
            *file_stream << HEADER_ENRICH_GENE_ID << FileSystem::DELIM_TSV
                         << HEADER_ENRICH_GO << std::endl;
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

/**
 * ======================================================================
 * Function bool QueryData::generate_transcriptome(uint32 flags, std::string outpath,
 *                                                  SEQUENCE_TYPES sequence_type)
 *
 * Description          - Generates a transcriptome with sequences based on the
 *                        flags input. If any flags match , it will be printed
 *
 * Notes                - None
 *
 * @param flags         - QuerySequence flags we would like to print
 * @param outpath       - Absolute outpath to print transcriptome
 * @param sequence_type - Sequence type we would like to print to the file
 *
 * @return              - Successful (true) or unsuccessful (false) generation of file
 *
 * =====================================================================
 */
bool QueryData::print_transcriptome(uint32 flags, std::string &outpath, SEQUENCE_TYPES sequence_type) {

    std::ofstream outfile(outpath, std::ios::out | std::ios::app);
    bool ret = true;
    uint64 sequence_ct = 0;

    try {
        // Ensure that we were able to open file for reading
        if (outfile.is_open()) {
            for (auto& pair : *mpSequences) {
                // If sequence flags match what we are looking for
                if (pair.second->getMQueryFlags() & flags) {
                    // Yes, print to the file
                    switch (sequence_type) {

                        case SEQUENCE_AMINO_ACID:
                            if (pair.second->is_protein()) {
                                outfile << pair.second->get_sequence_p() << std::endl;
                                sequence_ct++;
                            } else {
                                ;
                            }
                            break;

                        case SEQUENCE_NUCLEOTIDE:
                            if (pair.second->is_nucleotide()) {
                                outfile << pair.second->get_sequence_n() << std::endl;
                                sequence_ct++;
                            } else {
                                ;
                            }
                            break;

                        default:
                            FS_dprint("WARNING: QueryData print_transcriptome unrecognized case");
                            break;
                    }
                }
            }
        } else {
            // NO could not open file
            FS_dprint("ERROR unable to open file: " + outpath + " for writing");
            ret = false;
        };

    } catch (...) {
        ret = false;
    }

    if (sequence_ct == 0) {
        FS_dprint("WARNING: No sequences were printed to file " + outpath);
        ret = false;
    }
    outfile.close();
    return ret;
}
