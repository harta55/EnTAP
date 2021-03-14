/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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


#ifndef ENTAP_QUERYDATA_H
#define ENTAP_QUERYDATA_H


#include "QuerySequence.h"
#include "common.h"
#include "UserInput.h"

// Forward Declarations
class QueryAlignment;

typedef std::unordered_map<std::string, QuerySequence*> QUERY_MAP_T;

class QueryData {

public:

#ifdef UNIT_TESTS
    QueryData() = default;
#endif

    typedef enum {
        SUCCESS_EXPRESSION = (1 << 0),
        SUCCESS_FRAME_SEL  = (1 << 1),
        SUCCESS_ONTOLOGY   = (1 << 2),
        SUCCESS_SIM_SEARCH = (1 << 3),
        IS_PROTEIN         = (1 << 4),
        UNIPROT_MATCH      = (1 << 5),

        DATA_FLAGS_MAX     = (1 << 31)
    } DATA_FLAGS;

    typedef enum {
        SEQUENCE_AMINO_ACID=0,
        SEQUENCE_NUCLEOTIDE,
    } SEQUENCE_TYPES;


    QueryData(std::string &input_path, UserInput* userInput, FileSystem* fileSystem);
    QueryData(std::string&, std::string&, UserInput*, FileSystem*);
    ~QueryData();

    QUERY_MAP_T* get_sequences_ptr();

    std::pair<uint16, uint16> calculate_N_vals(std::vector<uint16>&,uint64);
    std::string trim_sequence_header(std::string&, std::string);
    void final_statistics(std::string& outpath);

    // Output routines
    bool start_alignment_files(const std::string &base_path, const std::vector<ENTAP_HEADERS> &headers, const vect_uint16_t &go_levels,
                                const std::vector<FileSystem::ENT_FILE_TYPES> &types);
    bool end_alignment_files(std::string &base_path);
    bool add_alignment_data(std::string &base_path, QuerySequence *querySequence, QueryAlignment *alignment);
    QuerySequence* get_sequence(std::string&);
    bool print_transcriptome(uint32 flags, std::string &outpath, SEQUENCE_TYPES sequence_type);

    QUERY_MAP_T get_specific_sequences(uint32 flags);

    // DATA_FLAG routines
    bool is_protein_data();
    void set_is_protein_data(bool val);
    void set_is_success_frame_selection(bool val);
    void set_is_success_expression(bool val);
    void set_is_success_sim_search(bool val);
    void set_is_success_ontology(bool val);
    void set_is_uniprot(bool val);

    // Header routines
    void header_set(ENTAP_HEADERS header, bool val);

#ifndef UNIT_TESTS
private:
#endif
    struct EntapHeader {
        const std::string title;
        bool print_header;
    };

    struct OutputFileData {
        std::vector<FileSystem::ENT_FILE_TYPES> file_types;
        vect_uint16_t go_levels;
        std::vector<ENTAP_HEADERS> headers;
        std::ofstream* file_streams[FileSystem::ENT_FILE_OUTPUT_FORMAT_MAX][UserInput::MAX_GO_LEVEL+1];
    };

    void set_input_type(std::string&);
    bool DATA_FLAG_GET(DATA_FLAGS);
    void DATA_FLAG_SET(DATA_FLAGS);
    void DATA_FLAG_CLEAR(DATA_FLAGS);
    void DATA_FLAG_CHANGE(DATA_FLAGS flag, bool val);

    void init_params(FileSystem *fileSystem, UserInput *userInput);
    std::string get_delim_data_sequence(std::vector<ENTAP_HEADERS>&headers, char delim, uint8 lvl,
                                QuerySequence* sequence);
    std::string get_delim_data_alignment(std::vector<ENTAP_HEADERS>&headers, char delim, uint8 lvl,
                                        QueryAlignment* alignment);
    bool initialize_file(std::ofstream *file_stream, std::vector<ENTAP_HEADERS> &headers, FileSystem::ENT_FILE_TYPES type);
    void generate_transcriptome(std::string &input_path, bool print_output, std::string output_path);

    const uint8         LINE_COUNT   = 20;
    const uint8         SEQ_DPRINT_CONUT = 10;
    const uint8         NUCLEO_DEV   = 2;
    const uint16        DEFAULT_GO_LEVEL = 0;           // Used if specified file does not have alternatives for go levels
    const fp32          N_50_PERCENT = 0.5;
    const fp32          N_90_PERCENT = 0.9;
    const std::string   NUCLEO_FLAG  = "Nucleotide";
    const std::string   PROTEIN_FLAG = "Protein";
    const std::string   COMPLETE_FLAG= "Complete";

    // these characters are considered nucleotide
    const std::vector<char> NUCLEO_MAP {
            'A', 'G', 'C', 'T', 'N'
    };
    const std::string OUT_UNANNOTATED_NUCL = "final_unannotated.fnn";
    const std::string OUT_UNANNOTATED_PROT = "final_unannotated.faa";
    const std::string OUT_ANNOTATED_NUCL   = "final_annotated.fnn";
    const std::string OUT_ANNOTATED_PROT   = "final_annotated.faa";
    const std::string APPEND_GO_LEVEL_STR  = "_lvl";
    const std::string APPEND_ENRICH_GENE_ID_GO = "_enrich_geneid_go";
    const std::string APPEND_ENRICH_GENE_ID_LEN = "_enrich_geneid_len";
    const std::string HEADER_ENRICH_GENE_ID = "gene_id";
    const std::string HEADER_ENRICH_GO      = "go_term";
    const std::string HEADER_ENRICH_LENGTH  = "effective_length";

    QUERY_MAP_T  *mpSequences;
    bool         mNoTrim;
    bool         mIsComplete;              // All sequences can be tagged as 'complete' genes
    uint32       mTotalSequences;          // Original sequence number
    uint32       mDataFlags;
    uint64       mNucleoLengthStart;       // Starting total len (nucleotide)
    uint64       mProteinLengthStart;      // Starting total len (protein)
    FileSystem  *mpFileSystem;
    UserInput   *mpUserInput;
    std::string mTranscriptTypeStr;
    std::unordered_map<std::string, OutputFileData> mAlignmentFiles;
    static EntapHeader ENTAP_HEADER_INFO[];

};


#endif //ENTAP_QUERYDATA_H
