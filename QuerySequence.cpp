//
// Created by harta on 3/29/17.
//

#include "QuerySequence.h"
#include <iostream>
#include <jmorecfg.h>

// best hit selection
bool QuerySequence::operator>(const QuerySequence &querySequence) {
    return this->e_val < querySequence.getE_val();

}

QuerySequence::QuerySequence(std::string database,std::string qseqid,std::string sseqid,
                             std::string stitle,float evalue) {
    this->database_path = database;
    this->qseqid = qseqid;
    this->sseqid = sseqid;
    this->stitle = stitle;
    this->e_val = evalue;
}

QuerySequence::QuerySequence() {

}

float QuerySequence::getE_val() const {
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
    return ostream << query.qseqid + "(" + query.sseqid + ", " << query.e_val << ", " + query.species + ")";
}
