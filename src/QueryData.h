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


#ifndef ENTAP_QUERYDATA_H
#define ENTAP_QUERYDATA_H


#include "QuerySequence.h"
#include "common.h"

// Forward Declarations
class QueryAlignment;

typedef std::unordered_map<std::string, QuerySequence*> QUERY_MAP_T;

class QueryData {

public:

    typedef enum {
        SUCCESS_EXPRESSION = (1 << 0),
        SUCCESS_FRAME_SEL  = (1 << 1),
        SUCCESS_ONTOLOGY   = (1 << 2),
        SUCCESS_SIM_SEARCH = (1 << 3),
        IS_PROTEIN         = (1 << 4),
        UNIPROT_MATCH      = (1 << 5),

        DATA_FLAGS_MAX     = (1 << 31)
    }DATA_FLAGS;


    QueryData(std::string&, std::string&, UserInput*, FileSystem*);
    ~QueryData();

    QUERY_MAP_T* get_sequences_ptr();

    std::pair<uint16, uint16> calculate_N_vals(std::vector<uint16>&,uint64);
    std::string trim_sequence_header(std::string&, std::string);
    void final_statistics(std::string& outpath);

    // Output routines
    bool start_alignment_files(std::string &base_path, std::vector<ENTAP_HEADERS> &headers, uint8 lvl,
                                std::vector<FileSystem::ENT_FILE_TYPES> &types);
    bool end_alignment_files(std::string &base_path);
    bool add_alignment_data(std::string &base_path, QuerySequence *querySequence, QueryAlignment *alignment);
    QuerySequence* get_sequence(std::string&);

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


private:

    struct EntapHeader {
        const std::string title;
        bool print_header;
    };

    struct OutputFileData {
        std::vector<FileSystem::ENT_FILE_TYPES> file_types;
        uint8 go_level;
        std::vector<ENTAP_HEADERS> headers;
        std::ofstream* file_streams[FileSystem::ENT_FILE_OUTPUT_FORMAT_MAX];
    };

    void set_input_type(std::string&);
    bool DATA_FLAG_GET(DATA_FLAGS);
    void DATA_FLAG_SET(DATA_FLAGS);
    void DATA_FLAG_CLEAR(DATA_FLAGS);
    void DATA_FLAG_CHANGE(DATA_FLAGS flag, bool val);

    std::string get_delim_data_sequence(std::vector<ENTAP_HEADERS>&headers, char delim, uint8 lvl,
                                QuerySequence* sequence);
    std::string get_delim_data_alignment(std::vector<ENTAP_HEADERS>&headers, char delim, uint8 lvl,
                                        QueryAlignment* alignment);
    bool initialize_file(std::ofstream *file_stream, std::vector<ENTAP_HEADERS> &headers, FileSystem::ENT_FILE_TYPES type);


    const uint8         LINE_COUNT   = 20;
    const uint8         SEQ_DPRINT_CONUT = 10;
    const uint8         NUCLEO_DEV   = 2;
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

    QUERY_MAP_T  *mpSequences;
    bool         mNoTrim;
    uint32       mTotalSequences;          // Original sequence number
    uint32       mDataFlags;
    uint64       mNucleoLengthStart;       // Starting total len (nucleotide)
    uint64       mProteinLengthStart;      // Starting total len (protein)
    uint32       mPipelineFlags;           // Success flags
    FileSystem  *mpFileSystem;
    UserInput   *mpUserInput;
    std::unordered_map<std::string, OutputFileData> mAlignmentFiles;
    static EntapHeader ENTAP_HEADER_INFO[];

};


#endif //ENTAP_QUERYDATA_H
