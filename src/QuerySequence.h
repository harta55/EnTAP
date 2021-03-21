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


#ifndef ENTAP_QUERYSEQUENCE_H
#define ENTAP_QUERYSEQUENCE_H

#include "common.h"
#include "database/EntapDatabase.h"

class QueryAlignment;

/**
 * ======================================================================
 * @class QuerySequence
 *
 * Description          - Handles sequence specific routines and data
 *                        management
 *                      - Allows setting of data for each of the supported
 *                        EnTAP modules/software (GeneMarkS-T, DIAMOND, etc)
 *                      - Organizes data corresponding to EnTAP headers into
 *                        mHeaderInfo for access during printing
 *                      - AlignmentData nested structure handles all alignment
 *                        information from each of the modules with alignments
 *
 * ======================================================================
 */
class QuerySequence {

public:

    //**************** Public Typedefs/Enums *******************
    typedef std::vector<QueryAlignment*> align_database_hits_t;

    typedef enum {

        QUERY_BLAST_HIT         = (1 << 0),         // Aligned against sequence during Similarity Search
        QUERY_EGGNOG_HIT        = (1 << 1),         // Aligned against sequence during EggNOG
        QUERY_EXPRESSION_KEPT   = (1 << 2),         // Was not removed during Expression Analysis
        QUERY_FRAME_KEPT        = (1 << 3),         // Was not removed during Frame Selection
        QUERY_FAMILY_ASSIGNED   = (1 << 4),         // Family was assigned during EggNOG
        QUERY_ONE_KEGG          = (1 << 5),         // Sequence contains at least one KEGG term from ANY process
        QUERY_ONE_GO            = (1 << 6),         // Sequence contains at least one GO term from ANY process
        QUERY_INFORMATIVE       = (1 << 7),         // Sequence is determined as informative from Sim Search
        QUERY_INTERPRO          = (1 << 8),         // Matched against InterPro databases
        QUERY_IS_PROTEIN        = (1 << 9),         // Sequence has a corresponding protein sequence
        QUERY_IS_NUCLEOTIDE     = (1 << 10),        // Sequence has a corresponding nucleotide sequence
        QUERY_BLASTED           = (1 << 11),        // Sequence was BLASTED (Sim Search was ran)
        QUERY_CONTAMINANT       = (1 << 12),        // Sequence is determined as a contaminant from Sim Search
        QUERY_FAMILY_ONE_KEGG   = (1 << 13),        // Sequence contains at least one KEGG from EggNOG process
        QUERY_FAMILY_ONE_GO     = (1 << 14),        // Sequence contains at least one GO from EggNOG process
        QUERY_ONT_INTERPRO_GO   = (1 << 15),        // Sequence contains at least one GO from InterPro process
        QUERY_ONT_INTERPRO_PATHWAY = (1 << 16),     // Sequence contains at least one KEGG from InterPro process
        QUERY_ONT_BUSCO         = (1 << 17),

        QUERY_MAX               = (1 << 31)

    } QUERY_FLAGS;
    //**********************************************************

    struct EggnogResults {
        std::string              member_ogs;        // 0A01R@biNOG,0V8CP@meNOG (ortholgous groups)
        std::string              seed_ortholog;     // 34740.HMEL017225-PA
        std::string              seed_evalue;       // Pulled from DIAMOND run
        std::string              seed_score;        // Pulled from DIAMOND run
        std::string              seed_coverage;     // Pulled from DIAMOND run
        std::string              predicted_gene;    // Most common predicted gene (pname)
        std::string              tax_scope_lvl_max; // virNOG[6]
        std::string              tax_scope;         // virNOG
        std::string              tax_scope_readable;// Ascomycota
        std::string              pname;             // All predicted gene names
        std::string              name;
        std::string              bigg;
        std::string              kegg;
        std::string              og_key;            // Used for indexing into older SQL database (if using)
        std::string              description;       // Used for older version
        std::string              protein_domains;
        fp64                     seed_eval_raw;     // Used for finding best hit
        go_format_t              parsed_go;         // All go terms found
    };

    struct InterProResults {
        std::string             e_value;
        std::string             database_desc_id;
        std::string             database_type;
        std::string             interpro_desc_id;
        std::string             pathways;
        fp64                    e_value_raw;
        go_format_t             parsed_go;
    };

    struct BuscoResults {
        std::string status; // Either Complete or Missing
        std::string busco_id;
        std::string length_str;
        uint64      length;
        fp64        score;
        std::string score_str;
    };

    struct SimSearchResults {
        std::string                       length;
        std::string                       mismatch;
        std::string                       gapopen;
        std::string                       qstart;
        std::string                       qend;
        std::string                       sstart;
        std::string                       send;
        std::string                       pident;
        std::string                       bit_score;
        std::string                       e_val;
        std::string                       coverage;
        std::string                       database_path;
        std::string                       qseqid;
        std::string                       sseqid;
        std::string                       stitle;
        std::string                       species;
        std::string                       contam_type;
        std::string                       lineage;
        std::string                       yes_no_contam; // just for convenience
        std::string                       yes_no_inform;
        fp32                              tax_score;     // taxonomic score, may be based on parent
        fp64                              e_val_raw;
        fp64                              coverage_raw;
        bool                              contaminant;
        bool                              is_informative;
        UniprotEntry                      uniprot_info;
    };

    /**
     * ======================================================================
     * @struct AlignmentData - nested
     *
     * Description          - Controls all alignment data from EnTAP modules
     *                        such as DIAMOND (multiple alignments/hits from
     *                        homology)
     *                      - Data is organized into maps of vectors (QueryAligments)
     *                        keyed to the database associated with the alignment
     *                      - First element in an alignment vector is the best hit
     *                        for that database
     *                      - All statuses in parent QuerySequence are based upon
     *                        the 'overall' best hit across all databases
     *                      - Best hit algorithm is implemented by the QueryAlignment
     *                        class
     *
     * ======================================================================
     */
    struct AlignmentData {
        typedef std::unordered_map<std::string,align_database_hits_t> ALIGNMENT_DATA_T;

        ALIGNMENT_DATA_T sim_search_data[SIM_SOFTWARE_COUNT];
        ALIGNMENT_DATA_T ontology_data[ONT_SOFTWARE_COUNT];
        QuerySequence* querySequence;

        QueryAlignment* overall_alignment[EXECUTION_MAX][ONT_SOFTWARE_COUNT]{};

        struct sort_descending_database {
            bool operator () (QueryAlignment* first, QueryAlignment* second);
        };

        AlignmentData(QuerySequence* sequence);
        ~AlignmentData();

        void set_best_alignment(QueryAlignment *new_alignment);
        void update_best_hit(QueryAlignment* new_alignment);
        bool hit_database(ExecuteStates state, uint16 software, std::string &database);
        align_database_hits_t* get_database_ptr(ExecuteStates, uint16, std::string&);
        QueryAlignment* get_best_align_ptr(ExecuteStates, uint16 software, std::string database);
        ALIGNMENT_DATA_T* get_software_ptr(ExecuteStates state, uint16 software);
    };



    /* Public Functions */
    QuerySequence();
    QuerySequence(bool, std::string, std::string);
    ~QuerySequence();
    void setFrame(const std::string &frame);
    uint64 get_sequence_length() const;
    const std::string &getFrame() const;
    const std::string &get_sequence_p() const;
    void set_sequence_p(std::string &seq);
    const std::string &get_sequence_n() const;
    void set_sequence_n(std::string &_sequence_n);
    const std::string &get_sequence() const;
    void set_fpkm(fp32 fpkm);
    const std::string &getMSequenceID() const;
    void setMTPM(fp64 mTPM);

    fp32 getMEffectiveLength() const;

    uint32 getMQueryFlags() const;

    void set_blasted();
    bool is_kept();
    bool is_protein();
    bool is_nucleotide();
    bool is_kept_expression();
    bool QUERY_FLAG_GET(QUERY_FLAGS flag);
    void QUERY_FLAG_CLEAR(QUERY_FLAGS flag);
    void QUERY_FLAG_CHANGE(QUERY_FLAGS flag, bool val);
    bool QUERY_FLAG_CONTAINS(uint32 flags);
    bool is_contaminant();
#ifdef EGGNOG_MAPPER
    void set_eggnog_results(const EggnogResults&);
#endif
    // Alignemnt accession routines
    void add_alignment(ExecuteStates state, uint16 software, EggnogResults &results, std::string& database);
    void add_alignment(ExecuteStates state, uint16 software, SimSearchResults &results, std::string& database,std::string lineage);
    void add_alignment(ExecuteStates state, uint16 software, InterProResults &results, std::string& database);
    void add_alignment(ExecuteStates state, uint16 software, BuscoResults &results, std::string& database);
    QuerySequence::align_database_hits_t* get_database_hits(std::string& database,ExecuteStates state, uint16 software);

    std::string format_go_info(go_format_t &go_list, uint8 lvl);

    // Returns recast alignment pointer
    template<class T>
    T *get_best_hit_alignment(ExecuteStates state, uint16 software, std::string database) {
        return static_cast<T*>(mAlignmentData->get_best_align_ptr(state, software, database));
    }

    // Checks whether an alignment was found against specific atabase
    bool hit_database(ExecuteStates state, uint16 software, std::string database);
    void update_query_flags(ExecuteStates state, uint16 software);
    void get_header_data(std::string& data, ENTAP_HEADERS header, uint8 lvl);
    void set_header_data();
    go_format_t get_go_terms();
    bool contains_go_level(int16 level);

    fp32 getMFrameScore() const;

    void setMFrameScore(fp32 mFrameScore);

    void setMEffectiveLength(fp32 mEffectiveLength);

private:
    //****************** Private Functions *********************
    void init_sequence();
    uint64 calc_seq_length(std::string &,bool);
    void trim_sequence(std::string& sequence);
    void QUERY_FLAG_SET(QUERY_FLAGS flag);
    //**********************************************************

    //**************** Private Const Variables *****************
    static constexpr uint8 NUCLEO_PER_AMINO = 3;        // Nucleotide/bp per Amino Acid
    //**********************************************************

    //****************** Private Variables *********************
    fp32                              mFPKM;            // FPKM value from Expression Filtering
    fp64                              mTPM;             // TPM value from Expression Filtering
    fp32                              mEffectiveLength; // Effective length from expression filtering
    uint32                            mQueryFlags;      // Status flags for the sequence
    std::string                       mSequenceID;      // Sequence ID
    uint64                            mSequenceLength;  // Sequence length (nucleotide bp)
    std::string                       mSequenceProtein; // Protein sequence
    std::string                       mSequenceNucleo;  // Nucleotide sequence
    std::string                       mFrameType;       // Frame type from Frame Selection
    fp32                              mFrameScore;      // Frame selection score
#ifdef EGGNOG_MAPPER
    EggnogResults                     mEggnogResults;   // EggNOG mapper results
#endif
    AlignmentData                     *mAlignmentData;  // Alignment information
    std::string                       mHeaderInfo[ENTAP_HEADER_COUNT];  // Header mappings
    //**********************************************************
};


#endif //ENTAP_QUERYSEQUENCE_H
