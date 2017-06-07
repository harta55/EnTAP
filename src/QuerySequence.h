//
// Created by harta on 3/29/17.
//

#ifndef ENTAP_QUERYSEQUENCE_H
#define ENTAP_QUERYSEQUENCE_H

#include <iostream>
#include <vector>
#include <string>
#include "Ontology.h"


class QuerySequence {
public:
    typedef std::map<std::string,std::vector<std::string>> go_struct;
    bool operator>(const QuerySequence& querySequence);
    void set_sim_search_results(std::string,std::string,std::string, double,int, int, int, int,int,
                  int, int, double, double, double, std::string);
    QuerySequence();
    QuerySequence(bool, std::string);
    friend void operator+(const QuerySequence &);
    void setQseqid(const std::string &qseqid);
    void setSequence(const std::string&);
    // TODO switch to map results
    void set_eggnog_results(std::string,std::string,std::string,std::string,std::string,
                    std::string,std::string,std::string);
    const std::string &get_contam_type() const;
    void set_contam_type(const std::string &_contam_type);
    void set_is_informative(bool _is_informative);
    void setIs_better_hit(bool is_better_hit);
    bool isContaminant() const;
    void setContaminant(bool contaminant);
    void set_tax_score(int _tax_score);
    void set_is_database_hit(bool _is_database_hit);
    void set_ontology_results(std::map<std::string,std::string>);
    std::string print_final_results(short,const std::vector<std::string>&);
    void set_lineage(const std::string &_lineage);
    void set_go_parsed(const go_struct &_go_parsed);
    const std::string &getSequence() const;
    void setSeq_length(unsigned long seq_length);
    void setFrame(const std::string &frame);
    void setSpecies(const std::string &species);
    unsigned long getSeq_length() const;
    const std::string &getFrame() const;
    bool isIs_protein() const;
    const std::string &get_species() const;
    bool is_informative() const;


private:
    friend std::ostream& operator<<(std::ostream& , const QuerySequence&);
    bool _contaminant, is_protein, is_better_hit, _is_informative, _is_database_hit;
    std::vector<std::string> _go_terms, _kegg_terms;
    int _tax_id,_length, _mismatch, _gapopen, _qstart, _qend, _sstart, _send, _tax_score;
    double _pident,_bit_score, _e_val, _coverage;
    unsigned long _seq_length;
    std::string _database_path, _qseqid,_sseqid, _stitle, _species, _sequence, _frame, _contam_type,
            _seed_ortho,_seed_eval,_seed_score,_predicted_gene,_tax_scope, _ogs,_go_str,_kegg_str,
            _lineage;
    go_struct _go_parsed;
    void init_sequence();
    std::map<std::string,std::string> _ontology_results;
    bool verify_frame(const std::string&,const std::string&);
};


#endif //ENTAP_QUERYSEQUENCE_H
