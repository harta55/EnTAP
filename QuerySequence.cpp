//
// Created by harta on 3/29/17.
//

#include "QuerySequence.h"
#include <iostream>
#include <string>
#include <algorithm>

// best hit selection
bool QuerySequence::operator>(const QuerySequence &querySequence) {
    if (this->is_better_hit) {
        // For hits of the same database "better hit"
        double eval1 = this->e_val, eval2 = querySequence.e_val;
        if (eval1 == 0) eval1 = 1E-120;
        if (eval2 == 0) eval2 = 1E-120;

        if (fabs(log10(eval1) - log10(eval2)) < 7) {
            if (this->contaminant && !querySequence.contaminant) return false;
            if (!this->contaminant && querySequence.contaminant) return true;

            double coverage_dif = fabs(this->_coverage - querySequence._coverage);
            if (coverage_dif < 4) {
                if (!this->_is_informative) return false;
                if (!querySequence._is_informative) return true;
                return this->length > querySequence.length;
            } else {
                return this->_coverage > querySequence._coverage;
            }
        } else {
            return eval1 > eval2;
        }
    }else {
        // For overall best hits between databases "best hit"
        if (this->contaminant && !querySequence.contaminant) return false;
        if (!this->contaminant && querySequence.contaminant) return true;

        double coverage_dif = fabs(this->_coverage - querySequence._coverage);
        if (coverage_dif < 7) {
            if (!this->_is_informative) return false;
            if (!querySequence._is_informative) return true;
            return this->length > querySequence.length;
        } else {
            return this->_coverage > querySequence._coverage;
        }
    }
}

void operator+(const QuerySequence &querySequence) {
//    this->database_path         = querySequence.database_path;
//    this->qseqid                = querySequence.qseqid;
//    this->sseqid                = querySequence.sseqid;

}

// TODO switch to set_sim_search
void QuerySequence::set_sim_search_results(std::string database,std::string qseqid,std::string sseqid,
                             double pident,int length, int mismatch, int gap, int qstart,
                             int qend, int sstart, int send, double evalue, double bit, double cover,
                             std::string title) {
    this->database_path = database;
    this->qseqid = qseqid;
    this->sseqid = sseqid;
    this->pident = pident;
    this->length = length;
    this->mismatch = mismatch;
    this->gapopen = gap;
    this->qstart = qstart;
    this->qend = qend;
    this->sstart = sstart;
    this->send = send;
    this->e_val = evalue;
    this->bit_score = bit;
    this->stitle = title;
    this->_coverage = cover;
}


unsigned long QuerySequence::getSeq_length() const {
    return seq_length;
}

QuerySequence::QuerySequence() {

}

void QuerySequence::setSequence(const std::string &seq) {
    this->_is_database_hit = false;
    this->is_protein = true;
    this->sequence = seq;
    if (!seq.empty() && seq[seq.length()-1] == '\n') {
        this->sequence.pop_back();
    }
}

bool QuerySequence::isIs_protein() const {
    return is_protein;
}

QuerySequence::QuerySequence(bool is_protein, std::string seq){
    this->_is_database_hit = false;
    this->is_protein = is_protein;
    std::string sub = seq.substr(seq.find("\n")+1);
    long line_chars = std::count(sub.begin(),sub.end(),'\n');
    unsigned long seq_len = sub.length() - line_chars;
    this->seq_length = seq_len;
    this->sequence = seq;
    if (!seq.empty() && seq[seq.length()-1] == '\n') {
        this->sequence.pop_back();
    }


}

void QuerySequence::setE_val(float e_val) {
    QuerySequence::e_val = e_val;
}

const std::string &QuerySequence::getDatabase_path() const {
    return database_path;
}

void QuerySequence::setDatabase_path(const std::string &database_path) {
    QuerySequence::database_path = database_path;
}

const std::string &QuerySequence::getQseqid() const {
    return qseqid;
}

void QuerySequence::setQseqid(const std::string &qseqid) {
    QuerySequence::qseqid = qseqid;
}

const std::string &QuerySequence::getSseqid() const {
    return sseqid;
}

void QuerySequence::setSseqid(const std::string &sseqid) {
    QuerySequence::sseqid = sseqid;
}

const std::string &QuerySequence::getStitle() const {
    return stitle;
}

void QuerySequence::setStitle(const std::string &stitle) {
    QuerySequence::stitle = stitle;
}

const std::string &QuerySequence::getSpecies() const {
    return species;
}

void QuerySequence::setSpecies(const std::string &species) {
    QuerySequence::species = species;
}

bool QuerySequence::isContaminant() const {
    return contaminant;
}

void QuerySequence::setContaminant(bool contaminant) {
    QuerySequence::contaminant = contaminant;
}

int QuerySequence::getTax_id() const {
    return tax_id;
}

void QuerySequence::setTax_id(int tax_id) {
    QuerySequence::tax_id = tax_id;
}

std::ostream& operator<<(std::ostream &ostream, const QuerySequence &query) {
    return ostream << query.qseqid<<'\t'<<query.sseqid<<'\t'<<query.pident<<'\t'<<
                    query.length<<'\t'<<query.mismatch<<'\t'<<query.mismatch<<'\t'<<
                    query.gapopen<<'\t'<<query.qstart<<'\t'<<query.qend<<'\t'<<
                    query.sstart<<'\t'<<query.send<<'\t'<<query.e_val<<'\t'<< query._coverage<<"\t"<<
                    query.stitle<<'\t'<<query.species<<'\t'<<query.database_path<<'\t'<<
                    query.frame;
}

const std::string &QuerySequence::getFrame() const {
    return frame;
}


double QuerySequence::getE_val() const {
    return e_val;
}

const std::string &QuerySequence::getSequence() const {
    return sequence;
}

void QuerySequence::setFrame(const std::string &frame) {
    QuerySequence::frame = frame;
}

void QuerySequence::setSeq_length(unsigned long seq_length) {
    QuerySequence::seq_length = seq_length;
}

// This is reserved for individual sim search file filtering
// Best hit for each database
void QuerySequence::setIs_better_hit(bool is_better_hit) {
    QuerySequence::is_better_hit = is_better_hit;
}

bool QuerySequence::is_is_informative() const {
    return _is_informative;
}

void QuerySequence::set_is_informative(bool _is_informative) {
    QuerySequence::_is_informative = _is_informative;
}

const std::string &QuerySequence::get_contam_type() const {
    return _contam_type;
}

void QuerySequence::set_contam_type(const std::string &_contam_type) {
    QuerySequence::_contam_type = _contam_type;
}

bool QuerySequence::is_is_database_hit() const {
    return _is_database_hit;
}

void QuerySequence::set_is_database_hit(bool _is_database_hit) {
    QuerySequence::_is_database_hit = _is_database_hit;
}
