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

        double getE_val() const;
        void setE_val(float e_val);
        const std::string &getDatabase_path() const;
        void setDatabase_path(const std::string &database_path);
        const std::string &getQseqid() const;
        void setQseqid(const std::string &qseqid);
        const std::string &getSseqid() const;
        void setSseqid(const std::string &sseqid);
        const std::string &getStitle() const;
        void setStitle(const std::string &stitle);
        void setSequence(const std::string&);
        void set_eggnog_results(std::string,std::string,std::string,std::string,std::string,
                        std::string,std::string,std::string);
        const std::string &get_contam_type() const;
        void set_contam_type(const std::string &_contam_type);
        bool is_is_informative() const;
        void set_is_informative(bool _is_informative);
        void setIs_better_hit(bool is_better_hit);
        bool isContaminant() const;
        void setContaminant(bool contaminant);
        int getTax_id() const;
        bool is_is_database_hit() const;
        void set_is_database_hit(bool _is_database_hit);
        void setTax_id(int tax_id);
        std::string print_eggnog();

    private:
        friend std::ostream& operator<<(std::ostream& , const QuerySequence&);
        bool contaminant, is_protein, is_better_hit, _is_informative, _is_database_hit;
        std::vector<std::string> _go_terms, _kegg_terms;
        int _tax_id,length, mismatch, gapopen, qstart, qend, sstart, send;
        double pident,bit_score, e_val, _coverage;
        unsigned long seq_length;
        std::string database_path, qseqid,sseqid, stitle, species, sequence, frame, _contam_type,
            _seed_ortho,_seed_eval,_seed_score,_predicted_gene,_tax_scope, _ogs,_go_str,_kegg_str;
        go_struct _go_parsed;
public:
    const go_struct &get_go_parsed() const;

    void set_go_parsed(const go_struct &_go_parsed);


public:
    void setSeq_length(unsigned long seq_length);

public:
    void setFrame(const std::string &frame);

public:
    const std::string &getSequence() const;


public:
    const std::string &getSpecies() const;

    void setSpecies(const std::string &species);

    unsigned long getSeq_length() const;

    const std::string &getFrame() const;

    bool isIs_protein() const;


};


#endif //ENTAP_QUERYSEQUENCE_H
