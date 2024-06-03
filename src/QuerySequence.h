/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2024, Alexander Hart, Dr. Jill Wegrzyn
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
        QUERY_EGGNOG_HIT        = (1 << 1),         // Aligned against sequence during EggNOG (either DIAMOND or mapper)
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
        QUERY_ONT_BUSCO         = (1 << 17),        // Sequence has BUSCO data
        QUERY_HGT_CANDIDATE     = (1 << 18),        // Sequence is an HGT candidate (not necessarily confirmed as HGT)
                                                    //  This means the sequence aligned with the correct number of donor/recipient
                                                    //  databases to be considered an HGT candidate
        QUERY_HGT_CONFIRMED     = (1 << 19),        // Sequence confirmed as HGT gene after performing neighbor analysis
        QUERY_HGT_BLASTED       = (1 << 20),        // Sequence hit against at least one donor or recipient database

        QUERY_MAX               = (1 << 31)

    } QUERY_FLAGS;
    //**********************************************************

    // With Eggnog, some data is shared between mapper and diamond implementations
    //  some are also not shared and only pulled from one version, outlined below
    struct EggnogResults {
        // Shared data
        std::string              seed_ortholog;     // 34740.HMEL017225-PA
        std::string              seed_evalue;       // Pulled from DIAMOND run
        std::string              seed_score;        // Pulled from DIAMOND run
        fp64                     seed_eval_raw;     // Used for finding best hit
        std::string              member_ogs;        // 0A01R@biNOG,0V8CP@meNOG (eggnog ortholgous groups)
        std::string              tax_scope_lvl_max; // 'virNOG[6]' or '33208|Metazoa' if using mapper
        std::string              description;       // Description of narrowest OG with a valid one
        go_format_t              parsed_go;         // All go terms found parsed into EnTAP format
        std::string              name;              // Preferred name
        std::string              bigg;              // BiGG reaction
        std::string              protein_domains;   // Pfam 'GCFC,NTR2' (comma separated)

        // DIAMOND specific eggnog data
        std::string              seed_coverage;     // Pulled from DIAMOND run
        std::string              predicted_gene;    // Most common predicted gene (pname)
        std::string              kegg;              // Everything combined in comma separated list
        std::string              tax_scope;         // virNOG
        std::string              tax_scope_readable;// Ascomycota
        std::string              pname;             // All predicted gene names
        std::string              og_key;            // Used for indexing into older SQL database (if using)

        // EggNOG-mapper specific eggnog data
        std::string              cog_category;      // COG category of narrowest OG with a valid one
        std::string              cog_category_description; // COG category description mapped from cog_category abbreviation
                                                           // EggNOG results may have multiple categories ('AJ' format)
                                                           // This will be formatted as 'description;description;'etc
        std::string              ec_value;
        std::string              kegg_ko;           // 'ko:K01672,ko:K12345' (comma separated)
        std::string              kegg_pathway;      // 'ko:K01672,ko:K12345' (comma separated)
        std::string              kegg_module ;      // 'M00355,M00595' (comma separated)
        std::string              kegg_reaction;     // 'R00132,R00154' (comma separated)
        std::string              kegg_rclass;       // 'RC02807,RC00299' (comma separated)
        std::string              brite;             // 'ko00000,ko03029 (comma separated)
        std::string              kegg_tc;           // '3.A.3.1,3.A.1.1' (comma separated)
        std::string              cazy;              // GH84
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

    struct HorizontalGeneTransferResults {
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
        std::string                       lineage;
        fp32                              tax_score;     // taxonomic score, may be based on parent
        fp64                              e_val_raw;
        fp64                              coverage_raw;
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
        ALIGNMENT_DATA_T horizontal_gene_data[HGT_SOFTWARE_COUNT];
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
    bool QUERY_FLAG_GET(QUERY_FLAGS flag) const;
    void QUERY_FLAG_CLEAR(QUERY_FLAGS flag);
    void QUERY_FLAG_CHANGE(QUERY_FLAGS flag, bool val);
    bool QUERY_FLAG_CONTAINS(uint32 flags);
    bool is_contaminant();

    // Alignemnt accession routines
    void add_alignment(ExecuteStates state, uint16 software, EggnogResults &results, std::string& database);
    void add_alignment(ExecuteStates state, uint16 software, SimSearchResults &results, std::string& database,std::string lineage);
    void add_alignment(ExecuteStates state, uint16 software, InterProResults &results, std::string& database);
    void add_alignment(ExecuteStates state, uint16 software, BuscoResults &results, std::string& database);
    void add_alignment(ExecuteStates state, uint16 software, HorizontalGeneTransferResults &results, std::string &database);
    QuerySequence::align_database_hits_t* get_database_hits(std::string& database,ExecuteStates state, uint16 software);

    std::string format_go_info(go_format_t &go_list, uint8 lvl);

    // Returns recast alignment pointer
    template<class T>
    T *get_best_hit_alignment(ExecuteStates state, uint16 software, std::string database) {
        return static_cast<T*>(mAlignmentData->get_best_align_ptr(state, software, database));
    }

    // Checks whether an alignment was found against specific database
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
    uint32                            mDonorDatabaseHitCt;
public:
    uint32 getMDonorDatabaseHitCt() const;
    void setMDonorDatabaseHitCt(uint32 mDonorDatabaseHitCt);
    uint32 getMRecipientDatabaseHitCt() const;
    void setMRecipientDatabaseHitCt(uint32 mRecipientDatabaseHitCt);
    const QuerySequence *getMpUpstreamSequence() const;
    void setMpUpstreamSequence(const QuerySequence *mpUpstreamSequence);
    const QuerySequence *getMpDownstreamSequence() const;
    void setMpDownstreamSequence(const QuerySequence *mpDownstreamSequence);
    
private:
    // Count of at least one alignment against donor database
    uint32                            mRecipientDatabaseHitCt; // Count of at least one alignment against recip database
    std::string                       mSequenceID;      // Sequence ID
    uint64                            mSequenceLength;  // Sequence length (nucleotide bp)
    std::string                       mSequenceProtein; // Protein sequence
    std::string                       mSequenceNucleo;  // Nucleotide sequence
    std::string                       mFrameType;       // Frame type from Frame Selection
    fp32                              mFrameScore;      // Frame selection score
    AlignmentData                     *mAlignmentData;  // Alignment information
    std::string                       mHeaderInfo[ENTAP_HEADER_COUNT];  // Header mappings

    /* Values taken from GFF file if user inputs */
    const QuerySequence *mpUpstreamSequence;
    const QuerySequence *mpDownstreamSequence;  // Sequence that is downstream from this sequence
};


#endif //ENTAP_QUERYSEQUENCE_H
