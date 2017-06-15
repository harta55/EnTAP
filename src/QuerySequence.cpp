//
// Created by harta on 3/29/17.
//

#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include "QuerySequence.h"
#include "EntapConsts.h"

// best hit selection
bool QuerySequence::operator>(const QuerySequence &querySequence) {
    if (this->is_better_hit) {
        // For hits of the same database "better hit"
        double eval1 = this->_e_val, eval2 = querySequence._e_val;
        if (eval1 == 0) eval1 = 1E-180;
        if (eval2 == 0) eval2 = 1E-180;
        if (fabs(log10(eval1) - log10(eval2)) < 6) {
            if (this->_contaminant && !querySequence._contaminant) return false;
            if (!this->_contaminant && querySequence._contaminant) return true;
            double coverage_dif = fabs(this->_coverage - querySequence._coverage);
            if (coverage_dif > 5) {
                return this->_coverage > querySequence._coverage;
            }
            if (this->_tax_score == querySequence._tax_score)
                return this->_e_val<querySequence._e_val;
            return this->_tax_score > querySequence._tax_score;
        } else {
            return eval1 < eval2;
        }
    }else {
        // For overall best hits between databases "best hit"
        if (this->_contaminant && !querySequence._contaminant) return false;
        if (!this->_contaminant && querySequence._contaminant) return true;
        double coverage_dif = fabs(this->_coverage - querySequence._coverage);
        if (coverage_dif > 5) {
            return this->_coverage > querySequence._coverage;
        }
        return this->_tax_score > querySequence._tax_score;
    }
}

void operator+(const QuerySequence &querySequence) {


}

// TODO switch to set_sim_search
void QuerySequence::set_sim_search_results(std::string database,std::string qseqid,std::string sseqid,
                             double pident,int length, int mismatch, int gap, int qstart,
                             int qend, int sstart, int send, double evalue, double bit, double cover,
                             std::string title) {
    _database_path = database;
    _qseqid = qseqid;
    _sseqid = sseqid;
    _pident = pident;
    _length = length;
    _mismatch = mismatch;
    _gapopen = gap;
    _qstart = qstart;
    _qend = qend;
    _sstart = sstart;
    _send = send;
    _e_val = evalue;
    _bit_score = bit;
    _stitle = title;
    _coverage = cover;
}


unsigned long QuerySequence::getSeq_length() const {
    return _seq_length;
}

QuerySequence::QuerySequence() {
    init_sequence();

}

void QuerySequence::setSequence( std::string &seq) {
    _is_database_hit = false;
    is_protein = true;
    _sequence = seq;
    _seq_length = calc_seq_length(seq,true);
    if (!seq.empty() && seq[seq.length()-1] == '\n') {
        this->_sequence.pop_back();
    }
}

bool QuerySequence::isIs_protein() const {
    return is_protein;
}

QuerySequence::QuerySequence(bool is_protein, std::string seq){
    init_sequence();
    this->_is_database_hit = false;
    this->is_protein = is_protein;
    _seq_length = calc_seq_length(seq,is_protein);
    this->_sequence = seq;
    if (!seq.empty() && seq[seq.length()-1] == '\n') {
        this->_sequence.pop_back();
    }
}

unsigned long QuerySequence::calc_seq_length(std::string &seq,bool protein) {
    std::string sub = seq.substr(seq.find("\n")+1);
    long line_chars = std::count(sub.begin(),sub.end(),'\n');
    unsigned long seq_len = sub.length() - line_chars;
    if (protein) seq_len *= 3;
    return seq_len;
}


void QuerySequence::setQseqid(const std::string &qseqid) {
    QuerySequence::_qseqid = qseqid;
}

void QuerySequence::setSpecies(const std::string &species) {
    QuerySequence::_species = species;
}

bool QuerySequence::isContaminant() const {
    return _contaminant;
}

void QuerySequence::setContaminant(bool contaminant) {
    QuerySequence::_contaminant = contaminant;
}

std::ostream& operator<<(std::ostream &ostream, const QuerySequence &query) {
    return ostream<<query._qseqid<<'\t' <<query._sseqid<<'\t'  <<query._pident<<'\t'<<
                    query._length<<'\t' <<query._mismatch<<'\t'<<
                    query._gapopen<<'\t'<<query._qstart<<'\t'  <<query._qend<<'\t'<<
                    query._sstart<<'\t' <<query._send<<'\t'    <<query._e_val<<'\t'<< query._coverage<<"\t"<<
                    query._stitle<<'\t' <<query._species<<'\t' <<query._database_path<<'\t'<<
                    query._frame;
}

const std::string &QuerySequence::getFrame() const {
    return _frame;
}

const std::string &QuerySequence::getSequence() const {
    return _sequence;
}

void QuerySequence::setFrame(const std::string &frame) {
    QuerySequence::_frame = frame;
}

void QuerySequence::setSeq_length(unsigned long seq_length) {
    QuerySequence::_seq_length = seq_length;
}

// This is reserved for individual sim search file filtering
// Best hit for each database
void QuerySequence::setIs_better_hit(bool is_better_hit) {
    QuerySequence::is_better_hit = is_better_hit;
}

void QuerySequence::set_is_informative(bool _is_informative) {
    QuerySequence::_is_informative = _is_informative;
}

const std::string &QuerySequence::get_species() const {
    return _species;
}

const std::string &QuerySequence::get_contam_type() const {
    return _contam_type;
}

bool QuerySequence::is_informative() const {
    return _is_informative;
}

void QuerySequence::set_contam_type(const std::string &_contam_type) {
    QuerySequence::_contam_type = _contam_type;
}

void QuerySequence::set_is_database_hit(bool _is_database_hit) {
    QuerySequence::_is_database_hit = _is_database_hit;
}

void QuerySequence::set_eggnog_results(std::string seed_o, std::string seed_o_eval,
                                       std::string seed_score, std::string predicted,
                                       std::string go_terms, std::string kegg,
                                       std::string annotation_tax, std::string ogs) {
    this->_go_str = go_terms;
    this->_kegg_str = kegg;
    this->_seed_ortho = seed_o;
    this->_seed_eval = seed_o_eval;
    this->_seed_score = seed_score;
    this->_predicted_gene = predicted;
    this->_tax_scope = annotation_tax;
    this->_ogs = ogs;
    std::stringstream ss(go_terms);
    std::string temp;
    if (!go_terms.empty()) {
        while (ss >> temp) {
            this->_go_terms.push_back(temp);
            if (ss.peek() == ',')
                ss.ignore();
        }
    }
    if (!kegg.empty()) {
        std::stringstream keggs(kegg);
        while (ss >> temp) {
            this->_kegg_terms.push_back(temp);
            if (ss.peek() == ',')
                ss.ignore();
        }
    }
}

std::string QuerySequence::print_final_results(short flag,const std::vector<std::string>&headers, short lvl) {
    std::stringstream stream;

    switch (flag) {
        case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
            stream << this->_qseqid       <<'\t' <<this->_sseqid    <<'\t'<<this->_pident        <<'\t'<<
                      this->_length       <<'\t' <<this->_mismatch  <<'\t'<<
                      this->_gapopen      <<'\t' <<this->_qstart    <<'\t'<<this->_qend          <<'\t'<<
                      this->_sstart       <<'\t' <<this->_send      <<'\t'<<this->_e_val         <<'\t'<<
                      this->_coverage     <<'\t' <<this->_stitle    <<'\t'<<this->_species       <<'\t'<<
                      this->_database_path<<'\t' <<this->_frame     <<'\t'<<this->_seed_ortho    <<'\t'<<
                      this->_seed_eval    <<'\t' <<this->_seed_score<<'\t'<<this->_predicted_gene<<'\t'<<
                      this->_tax_scope    <<'\t'<< this->_ogs       <<'\t'<<this->_kegg_str      <<'\t';
            break;
        case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
            stream << *this <<'\t';
            for (const std::string &val : headers) {
                stream << _ontology_results[val] << '\t';
            }
            break;
        default:
            break;
    }
    if (!this->_go_str.empty()) {
        for (std::string val : _go_parsed[ENTAP_EXECUTE::GO_BIOLOGICAL_FLAG]) {
            if (val.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                stream<<val<<",";
            }
        }
        stream<<'\t';
        for (std::string val : _go_parsed[ENTAP_EXECUTE::GO_CELLULAR_FLAG])  {
            if (val.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                stream<<val<<",";
            }
        }
        stream<<'\t';
        for (std::string val : _go_parsed[ENTAP_EXECUTE::GO_MOLECULAR_FLAG]) {
            if (val.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                stream<<val<<",";
            }
        }
        stream<<'\t';
    } else {
        stream<<"\t\t\t";
    }
    return stream.str();

}

void QuerySequence::set_go_parsed(const QuerySequence::go_struct &_go_parsed) {
    QuerySequence::_go_parsed = _go_parsed;
}

void QuerySequence::init_sequence() {
    _seq_length = 0;
    _pident = 0;
    _length = 0;
    _mismatch = 0;
    _gapopen = 0;
    _qstart = 0;
    _qend = 0;
    _sstart = 0;
    _send = 0;
    _e_val = 0;
    _bit_score = 0;
    _coverage = 0;
}

void QuerySequence::set_ontology_results(std::map<std::string, std::string> map) {
    this->_ontology_results = map;
}

bool QuerySequence::verify_frame(const std::string &frame1, const std::string &frame2) {
    return (frame1.compare(ENTAP_EXECUTE::FRAME_SELECTION_INTERNAL_FLAG)!=0 &&
            frame2.compare(ENTAP_EXECUTE::FRAME_SELECTION_INTERNAL_FLAG)!=0);
}

void QuerySequence::set_lineage(const std::string &_lineage) {
    QuerySequence::_lineage = _lineage;
}

void QuerySequence::set_tax_score(int _tax_score) {
    QuerySequence::_tax_score = _tax_score;
};
