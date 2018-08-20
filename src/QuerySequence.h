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


#ifndef ENTAP_QUERYSEQUENCE_H
#define ENTAP_QUERYSEQUENCE_H

#include "common.h"
#include "Ontology.h"
#include "database/SQLDatabaseHelper.h"
#include "EntapExecute.h"
class QuerySequence {
public:

    class QueryAlignment;
    struct AlignmentData;

    // Avoid cluttering global with the following defs
    typedef std::pair<QueryAlignment*,std::vector<QueryAlignment*>> align_database_hits_t;
    typedef std::unordered_map<std::string,align_database_hits_t> ALIGNMENT_DATA_T;

    typedef enum {

        QUERY_BLAST_HIT         = (1 << 0),
        QUERY_EGGNOG_HIT        = (1 << 1),
        QUERY_EXPRESSION_KEPT   = (1 << 2),
        QUERY_FRAME_KEPT        = (1 << 3),
        QUERY_FAMILY_ASSIGNED   = (1 << 4),
        QUERY_ONE_KEGG          = (1 << 5),
        QUERY_ONE_GO            = (1 << 6),
        QUERY_INFORMATIVE       = (1 << 7),
        QUERY_INTERPRO          = (1 << 8),
        QUERY_IS_PROTEIN        = (1 << 9),
        QUERY_BLASTED           = (1 << 10),
        QUERY_CONTAMINANT       = (1 << 11)

    } QUERY_FLAGS;

    struct EggnogResults {
        std::string              best_hit_query; // 34740.HMEL017225-PA
        std::string              member_ogs;     // 0A01R@biNOG,0V8CP@meNOG
        std::string              seed_ortholog;
        std::string              seed_evalue;
        std::string              seed_score;
        std::string              predicted_gene;
        std::string              tax_scope_lvl_max; // virNOG[6]
        std::string              tax_scope;         // virNOG NOT virNOG[6]
        std::string              tax_scope_readable;// Ascomycota
        std::string              pname;
        std::string              name;
        std::string              bigg;
        std::string              ogs;
        std::string              og_key;
        std::string              sql_kegg;
        std::string              sql_go;
        std::string              description;
        std::string              protein_domains;
        std::vector<std::string> raw_kegg;
        std::vector<std::string> raw_go;
        go_format_t              parsed_go;
    };

    struct InterProResults {
        std::string             e_value;
        std::string             database_desc_id;
        std::string             database_type;
        std::string             interpro_desc_id;
        std::string             pathways;
        go_format_t             parsed_go;
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

    QuerySequence();
    QuerySequence(bool, std::string, std::string);
    ~QuerySequence();
    void setSequence(std::string&);
    std::string print_tsv(const std::vector<const std::string*>&);
    std::string print_tsv(std::vector<const std::string*>& , short);
    void init_header();
    void setSeq_length(unsigned long seq_length);
    void setFrame(const std::string &frame);
    unsigned long getSeq_length() const;
    const std::string &getFrame() const;
    const std::string &get_sequence_p() const;
    void set_sequence_p(const std::string &_sequence_p);
    const std::string &get_sequence_n() const;
    void set_sequence_n(const std::string &_sequence_n);
    const std::string &get_sequence() const;
    void set_fpkm(float _fpkm);
    bool is_kept();
    bool QUERY_FLAG_GET(QUERY_FLAGS);
    void QUERY_FLAG_SET(QUERY_FLAGS);
    void QUERY_FLAG_CLEAR(QUERY_FLAGS);
    bool isContaminant();
    bool hit_database(ExecuteStates state, uint16 software, std::string database);
    void set_eggnog_results(const EggnogResults&);
    void set_interpro_results(std::string&,std::string&,std::string&,std::string&,
                              std::string&,go_format_t&);
    QuerySequence::align_database_hits_t* get_database_hits(std::string&,ExecuteStates, uint16 );
    std::string get_header_data(const std::string*);

    template<class T>
    T *get_best_hit_alignment(ExecuteStates state, uint16 software, std::string database) {
        if (database.empty()) {
            return static_cast<T*>(_total_alignment_data.index_best_align(state,software));
        } else {
            return static_cast<T*>(_total_alignment_data.index_data(state,software)->at(database).first);
        }
    }

    template<class U>
    void add_alignment(ExecuteStates state, uint16 software, U results, std::string& database, std::string& lineage) {
        // Create new alignment object
        // Update vector containing all alignments
        QueryAlignment *new_alignment;
        switch (state) {
            case SIMILARITY_SEARCH: {
                QUERY_FLAG_SET(QUERY_BLAST_HIT);
                SimSearchResults cast_results = static_cast<SimSearchResults>(results);
                new_alignment = new SimSearchAlignment(cast_results, lineage, this);
                break;
            }
            case EXPRESSION_FILTERING:
            case FRAME_SELECTION:
            case FILTER:
            case GENE_ONTOLOGY:
            case EXIT:
            case EXECUTION_MAX:
            default:
                return;
        }

        // Did we not hit against the database yet?
        if (!hit_database(state, software, database)) {
            // No, create new vector for that database and add as best hit for database
            std::vector<QueryAlignment*> vect = {new_alignment};
            _total_alignment_data.index_data(state,software)->emplace(database, std::make_pair(new_alignment, vect));
        } else {
            // Yes, add alignment to list then update
            _total_alignment_data.index_data(state,software)->at(database).second.push_back(new_alignment);
        }

        // Update list and best hits
        update_best_hit(state, software, database);
    }

    void update_best_hit(ExecuteStates state, uint16 software, std::string &database);
public:

    //**********************************************************************
    //**********************************************************************
    //                 QueryAlignment Nested Class
    //**********************************************************************
    //**********************************************************************

    class QueryAlignment {

    public:
        QueryAlignment();
        std::string print_tsv(const std::vector<const std::string*>&);
        bool operator<(const QueryAlignment&query) {return !(*this > query);};
        void set_compare_overall_alignment(bool val);
        virtual ~QueryAlignment(){};
        virtual bool operator>(const QueryAlignment&)=0;

    protected:
        std::map<const std::string*, std::string> ALIGN_OUTPUT_MAP;
        bool _compare_overall_alignment; // May want to compare separate parameters for overall alignment
        QuerySequence* _parent;
    };

    class SimSearchAlignment : public QueryAlignment{

    public:
        SimSearchAlignment(SimSearchResults, std::string&, QuerySequence*);
        virtual ~SimSearchAlignment();
        SimSearchResults* get_results();
        bool operator>(const QueryAlignment&);

    private:
        void set_tax_score(std::string&);

        QuerySequence::SimSearchResults    _sim_search_results;

    protected:
        static constexpr uint8 E_VAL_DIF     = 8;
        static constexpr uint8 COV_DIF       = 5;
        static constexpr uint8 INFORM_ADD    = 3;
        static constexpr fp32 INFORM_FACTOR  = 1.2;
    };


    struct AlignmentData {
        std::vector<ALIGNMENT_DATA_T>     _alignment_data;  // contains all alignment data
        std::vector<QueryAlignment*>      _alignment_best;  // only overall best hits

        AlignmentData();
        ~AlignmentData();

        ALIGNMENT_DATA_T* index_data(ExecuteStates state, uint16 software);
        QueryAlignment* index_best_align(ExecuteStates state, uint16 software);
        void set_best_align(ExecuteStates state, uint16 software, QueryAlignment*);
    };

private:
    fp32                              _fpkm;
    uint16                            _query_flags;
    std::string                       _seq_id;
    unsigned long                     _seq_length;
    std::string                       _sequence_p;
    std::string                       _sequence_n;
    std::string                       _frame;
    EggnogResults                     _eggnog_results;
    InterProResults                   _interpro_results;
    std::map<const std::string*, std::string*> OUTPUT_MAP;
    AlignmentData                     _total_alignment_data;

    void init_sequence();
    unsigned long calc_seq_length(std::string &,bool);
    void update_query_flags(ExecuteStates state, uint16 software);
};


#endif //ENTAP_QUERYSEQUENCE_H
