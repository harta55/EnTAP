//
// Created by harta on 3/29/17.
//

#include "QuerySequence.h"
#include <iostream>

// best hit selection
bool QuerySequence::operator>(const QuerySequence &querySequence) {
    if (this->contaminant && !querySequence.contaminant) return false;
    if (!this->contaminant && querySequence.contaminant) return true;
    if (this->e_val<user_e && querySequence.e_val>user_e) return true;
    if (this->e_val>user_e && querySequence.e_val<user_e) return false;
    return this->length>querySequence.length;
}

QuerySequence::QuerySequence(std::string database,std::string qseqid,std::string sseqid,
                             float pident,int length, int mismatch, int gap, int qstart,
                             int qend, int sstart, int send, double evalue, float bit,
                             std::string title, double user_e) {
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
    this->user_e = user_e;
}

QuerySequence::QuerySequence() {

}

double QuerySequence::getE_val() const {
    return e_val;
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
                    query.sstart<<'\t'<<query.send<<'\t'<<query.e_val<<'\t'<<
                    query.informative<<'\t'<<query.species<<'\t'<<query.database_path;
}

const std::string &QuerySequence::getInformative() const {
    return informative;
}

void QuerySequence::setInformative(const std::string &informative) {
    QuerySequence::informative = informative;
}
