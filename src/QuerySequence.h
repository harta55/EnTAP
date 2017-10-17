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


#ifndef ENTAP_QUERYSEQUENCE_H
#define ENTAP_QUERYSEQUENCE_H

#include <iostream>
#include <vector>
#include <string>
#include "Ontology.h"
#include "DatabaseHelper.h"

class QuerySequence {
public:

    typedef std::map<std::string,std::vector<std::string>> go_struct;

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
        go_struct                parsed_go;
    };

    struct InterProResults {
        std::string             e_value;
        std::string             database_desc_id;
        std::string             database_type;
        std::string             interpro_desc_id;
        std::string             pathways;
        go_struct               parsed_go;
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
        bool                              contaminant;
    };

    bool operator>(const QuerySequence& querySequence);
    void set_sim_search_results(std::string,std::string,std::string,
                                std::string,std::string, std::string, std::string, std::string,
                                std::string, std::string, std::string, std::string, std::string,
                                double,  double);
    QuerySequence();
    QuerySequence(bool, std::string);
    friend void operator+(const QuerySequence &);
    void setQseqid(const std::string &qseqid);
    void setSequence(std::string&);
    // TODO switch to map results
    void set_eggnog_results(std::string,std::string,std::string,std::string,std::string,
                    std::string,std::string,std::string, DatabaseHelper &);
    void set_interpro_results(std::string&,std::string&,std::string&,std::string&,
                              std::string&,go_struct&);
    std::string print_tsv(const std::vector<const std::string*>&);
    std::string print_tsv(std::vector<const std::string*>& , short);
    std::string format_eggnog(std::string&);

    void set_tax_score(std::string);
    void init_header();
    const std::string &get_contam_type() const;
    const SimSearchResults &get_sim_struct() const;
    void set_sim_struct(const SimSearchResults &);
    void set_contam_type(const std::string &_contam_type);
    void set_is_informative(bool _is_informative);
    void setIs_better_hit(bool is_better_hit);
    bool isContaminant() const;
    void setContaminant(bool contaminant);
    void set_is_database_hit(bool _is_database_hit);
    void set_lineage(const std::string &_lineage);
    void set_go_parsed(const go_struct &_go_parsed, short);
    void setSeq_length(unsigned long seq_length);
    void setFrame(const std::string &frame);
    void setSpecies(const std::string &species);
    unsigned long getSeq_length() const;
    const std::string &getFrame() const;
    bool isIs_protein() const;
    const std::string &get_species() const;
    bool is_informative() const;
    const std::string &get_sequence_p() const;
    void set_sequence_p(const std::string &_sequence_p);
    const std::string &get_sequence_n() const;
    void set_sequence_n(const std::string &_sequence_n);
    const std::string &get_sequence() const;
    void setIs_protein(bool is_protein);
    bool is_is_database_hit() const;
    bool is_is_family_assigned() const;
    void set_is_family_assigned(bool _is_family_assigned);
    bool is_is_one_go() const;
    void set_is_one_go(bool _is_one_go);
    bool is_is_one_kegg() const;
    void set_is_one_kegg(bool _is_one_kegg);
    bool is_is_expression_kept() const;
    void set_is_expression_kept(bool _is_expression_kept);
    void set_fpkm(float _fpkm);
    const std::string &get_tax_scope() const;
    void set_is_interpro_hit(bool _is_interpro_hit);
    bool is_kept() const;
    void set_kept(bool _kept);

private:

    static constexpr uint8 E_VAL_DIF     = 8;
    static constexpr uint8 COV_DIF       = 5;
    static constexpr uint8 INFORM_ADD    = 3;
    static constexpr fp32 INFORM_FACTOR  = 1.2;

    bool                              _kept;        // Continually updated
    bool                              is_protein;
    bool                              is_better_hit;
    bool                              _is_informative;
    bool                              _is_database_hit;
    bool                              _is_family_assigned;
    bool                              _is_one_go;
    bool                              _is_one_kegg;
    bool                              _is_expression_kept;
    bool                              _is_interpro_hit;
    bool                              _is_eggnog_hit;
    float                             _tax_score;
    float                             _fpkm;
    unsigned long                     _seq_length;
    double                            _e_val;
    double                            _coverage;
    std::string                       _sequence_p;
    std::string                       _sequence_n;
    std::string                       _frame;
    EggnogResults                     _eggnog_results;
    SimSearchResults                  _sim_search_results;
    InterProResults                   _interpro_results;
    std::map<const std::string*, std::string*> OUTPUT_MAP;

    friend std::ostream& operator<<(std::ostream& , const QuerySequence&);
    void init_sequence();
    unsigned long calc_seq_length(std::string &,bool);
};


#endif //ENTAP_QUERYSEQUENCE_H
