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

#include <iostream>
#include <vector>
#include <string>
#include "Ontology.h"
#include "database/SQLDatabaseHelper.h"
#include "EntapExecute.h"


class QuerySequence;
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
    };
// Taken from best hit

class QueryAlignment {

public:
    QueryAlignment(QuerySequence*, uint16, std::string);
    virtual std::string print_tsv(const std::vector<const std::string*>&)=0;
    virtual ~QueryAlignment() = default;

protected:
    QuerySequence   *_parent_sequence;
    uint16           _software_flag;
    std::string      _database;
    std::string      _frame;
    std::map<const std::string*, std::string*> OUTPUT_MAP;
};

class SimSearchAlignment : public QueryAlignment{

public:
    bool operator>(const SimSearchAlignment&);
    bool operator<(const SimSearchAlignment&query) {return !(*this > query);};
    SimSearchAlignment(QuerySequence*, uint16, std::string, SimSearchResults,
                       std::string&);
    SimSearchResults* get_results();
    std::string print_tsv(const std::vector<const std::string*>&header){ return QueryAlignment::print_tsv(header); }
    void set_best_hit(bool a){this->_best_hit = a;}

private:
    void set_tax_score(std::string&);

    SimSearchResults    _sim_search_results;
    bool                _best_hit;

protected:
    static constexpr uint8 E_VAL_DIF     = 8;
    static constexpr uint8 COV_DIF       = 5;
    static constexpr uint8 INFORM_ADD    = 3;
    static constexpr fp32 INFORM_FACTOR  = 1.2;
};

class QuerySequence {
public:

    // Avoid cluttering global with the following defs
    typedef std::pair<QueryAlignment*,std::vector<QueryAlignment*>> align_database_hits_t;
    typedef std::unordered_map<std::string,align_database_hits_t> alignment_database_map_t;

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
        std::string              seed_ortholog;
        std::string              seed_evalue;
        std::string              seed_score;
        std::string              predicted_gene;
        std::string              tax_scope;         // virNOG NOT virNOG[6]
        std::string              tax_scope_readable;// Ascomycota
        std::string              ogs;
        std::string              og_key;
        std::string              sql_kegg;
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

    struct SimSearchAlignmentData {
        alignment_database_map_t alignments;
        SimSearchAlignment* best_hit;
        SimSearchResults *results;

        SimSearchAlignmentData() {
            alignments = {};
            best_hit = nullptr;
            results = nullptr;
        }
    };

    QuerySequence();
    QuerySequence(bool, std::string, std::string);
    ~QuerySequence();
    void setSequence(std::string&);
    // TODO switch to map results
    void set_eggnog_results(const EggnogResults&);
    void set_interpro_results(std::string&,std::string&,std::string&,std::string&,
                              std::string&,go_format_t&);
    std::string print_tsv(const std::vector<const std::string*>&);
    std::string print_tsv(std::vector<const std::string*>& , short);

    void init_header();
    const std::string &get_contam_type() const;
    void setSeq_length(unsigned long seq_length);
    void setFrame(const std::string &frame);

    unsigned long getSeq_length() const;
    const std::string &getFrame() const;
    const std::string &get_species() const;
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

    bool hit_database(std::string&, ExecuteStates);


    template<class T>
    T *get_best_hit_alignment(ExecuteStates state, std::string database) {
        switch (state) {
            case SIMILARITY_SEARCH:
                if (database.empty()) {
                    return static_cast<T*>(this->_sim_search_alignment_data->best_hit);
                }else {
                    return static_cast<T*>(this->_sim_search_alignment_data->alignments[database].first);
                }
            default:
                return nullptr;
        }
    }

    template<class T>
    T *set_best_hit_alignment(ExecuteStates state, std::string database, QueryAlignment* alignment) {

        switch (state) {
            case SIMILARITY_SEARCH:
                if (database.empty()) {
                    this->_sim_search_alignment_data->best_hit = static_cast<T*>(alignment);
                }else {
                    this->_sim_search_alignment_data->alignments[database].first = static_cast<T*>(alignment);
                }
            default:
                return nullptr;
        }
    }

    template<class T, class U>
    void add_alignment(ExecuteStates state, uint16 software, U &results, std::string &database,
                                      std::string &lineage) {
        // Create new alignment object
        T *new_alignment = (new T(this, software, database, results, lineage));
        // Update vector containing all alignments
        switch (state) {
            case SIMILARITY_SEARCH:
                QUERY_FLAG_SET(QUERY_BLAST_HIT);
                if (this->_sim_search_alignment_data->alignments.find(database) ==
                    this->_sim_search_alignment_data->alignments.end()) {
                    // database not found, make new entry for database
                    std::vector<QueryAlignment*> vect = {new_alignment};
                    this->_sim_search_alignment_data->alignments.emplace(database,
                            std::make_pair(nullptr, vect));
                }else {
                    // database found, just add to vector containing all hits
                    this->_sim_search_alignment_data->alignments[database].second.push_back(
                            new_alignment
                    );
                }
                break;
            default:
                return;
        }
        update_best_hit<T>(state, database);
    }

    template<class T>
    void update_best_hit(ExecuteStates state, std::string &database) {
        align_database_hits_t *database_data = get_database_hits(database, state);
        // first: database, second: pair of best hit to all hits (vector)
        for (auto *temp : database_data->second) {
            T* alignment = static_cast<T*>(temp);
            T* database_best_hit = get_best_hit_alignment<T>(state, database);
            // cycle through alignments and update the best hit
            if (database_best_hit == nullptr || *alignment > *database_best_hit) {
                // Update best hit for database
                set_best_hit_alignment<T>(state,database,temp);
                // Update overall best hit for parent
                T* best_overall_hit = get_best_hit_alignment<T>(state, "");
                switch (state) {
                    case SIMILARITY_SEARCH:
                        // must only compare best overall (between databases) for sim search)=
                        alignment->set_best_hit(true);
                        break;
                    default:
                        break;
                }
                if (best_overall_hit == nullptr || *alignment > *best_overall_hit) {
                    set_best_hit_alignment<T>(state, "", temp);
                    // Update parent flags, data...
                    update_query_flags(state);
                }
                switch (state) {
                    case SIMILARITY_SEARCH:
                        // must only compare best overall (between databases) for sim search)=
                        alignment->set_best_hit(false);
                        break;
                    default:
                        break;
                }
            }
        }
    }

    QuerySequence::align_database_hits_t* get_database_hits(std::string&, ExecuteStates);

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
    SimSearchAlignmentData            *_sim_search_alignment_data;  // contains all alignment data
    std::map<const std::string*, std::string*> OUTPUT_MAP;

    void init_sequence();
    unsigned long calc_seq_length(std::string &,bool);
    void update_query_flags(ExecuteStates);
};


#endif //ENTAP_QUERYSEQUENCE_H
