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


#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include "QuerySequence.h"
#include "EntapGlobals.h"
#include "EggnogLevels.h"

// best hit selection
bool QuerySequence::operator>(const QuerySequence &querySequence) {
    if (this->is_better_hit) {
        // For hits of the same database "better hit"
        double eval1 = this->_e_val, eval2 = querySequence._e_val;
        if (eval1 == 0) eval1 = 1E-180;
        if (eval2 == 0) eval2 = 1E-180;
        if (fabs(log10(eval1) - log10(eval2)) < E_VAL_DIF) {
            double coverage_dif = fabs(this->_coverage - querySequence._coverage);
            if (coverage_dif > COV_DIF) {
                return this->_coverage > querySequence._coverage;
            }
            if (this->_sim_search_results.contaminant && !querySequence._sim_search_results.contaminant) return false;
            if (!this->_sim_search_results.contaminant && querySequence._sim_search_results.contaminant) return true;
            if (this->_tax_score == querySequence._tax_score)
                return this->_e_val<querySequence._e_val;
            return this->_tax_score > querySequence._tax_score;
        } else {
            return eval1 < eval2;
        }
    }else {
        // For overall best hits between databases "best hit"
        double coverage_dif = fabs(this->_coverage - querySequence._coverage);
        if (coverage_dif > COV_DIF) {
            return this->_coverage > querySequence._coverage;
        }
        if (this->_sim_search_results.contaminant && !querySequence._sim_search_results.contaminant) return false;
        if (!this->_sim_search_results.contaminant && querySequence._sim_search_results.contaminant) return true;
        return this->_tax_score > querySequence._tax_score;
    }
}

void operator+(const QuerySequence &querySequence) {


}

// TODO switch to set_sim_search
void QuerySequence::set_sim_search_results(std::string database,std::string qseqid,std::string sseqid,
                             std::string pident,std::string length, std::string mismatch, std::string gap, std::string qstart,
                             std::string qend , std::string sstart, std::string send, std::string title, std::string bit,
                             double evalue,  double cover) {
    _sim_search_results.database_path = database;
    _sim_search_results.qseqid = qseqid;
    _sim_search_results.sseqid = sseqid;
    _sim_search_results.pident = pident;
    _sim_search_results.length = length;
    _sim_search_results.mismatch = mismatch;
    _sim_search_results.gapopen = gap;
    _sim_search_results.qstart = qstart;
    _sim_search_results.qend = qend;
    _sim_search_results.sstart = sstart;
    _sim_search_results.send = send;
    _sim_search_results.stitle = title;
    _sim_search_results.bit_score = bit;
    std::ostringstream ostringstream;
    ostringstream<<evalue;
    _sim_search_results.e_val = ostringstream.str();
    _sim_search_results.coverage = std::to_string(cover);

    _e_val = evalue;
    _coverage = cover;
}


unsigned long QuerySequence::getSeq_length() const {
    return _seq_length;
}

QuerySequence::QuerySequence() {
    init_sequence();

}

void QuerySequence::setSequence( std::string &seq) {
    is_protein = true;
    _seq_length = calc_seq_length(seq,true);
    if (!seq.empty() && seq[seq.length()-1] == '\n') {
        seq.pop_back();
    }
    _sequence_p = seq;
}

const std::string &QuerySequence::get_sequence_p() const {
    return _sequence_p;
}

void QuerySequence::set_sequence_p(const std::string &_sequence_p) {
    QuerySequence::_sequence_p = _sequence_p;
}

const std::string &QuerySequence::get_sequence_n() const {
    return _sequence_n;
}

void QuerySequence::set_sequence_n(const std::string &_sequence_n) {
    QuerySequence::_sequence_n = _sequence_n;
}

bool QuerySequence::isIs_protein() const {
    return is_protein;
}

QuerySequence::QuerySequence(bool is_protein, std::string seq){
    init_sequence();
    this->_is_database_hit = false;
    this->is_protein = is_protein;
    _seq_length = calc_seq_length(seq,is_protein);
    if (!seq.empty() && seq[seq.length()-1] == '\n') {
        seq.pop_back();
    }
    is_protein ? _sequence_p = seq : _sequence_n = seq;
}

unsigned long QuerySequence::calc_seq_length(std::string &seq,bool protein) {
    std::string sub = seq.substr(seq.find("\n")+1);
    long line_chars = std::count(sub.begin(),sub.end(),'\n');
    unsigned long seq_len = sub.length() - line_chars;
    if (protein) seq_len *= 3;
    return seq_len;
}


void QuerySequence::setQseqid(const std::string &qseqid) {
    _sim_search_results.qseqid = qseqid;
}

void QuerySequence::setSpecies(const std::string &species) {
    _sim_search_results.species = species;
}

bool QuerySequence::isContaminant() const {
    return _sim_search_results.contaminant;
}

void QuerySequence::setContaminant(bool contaminant) {
    QuerySequence::_sim_search_results.contaminant = contaminant;
    contaminant ? _sim_search_results.yes_no_contam = "Yes" : _sim_search_results.yes_no_contam  = "No";
}

std::ostream& operator<<(std::ostream &ostream, const QuerySequence &query) {
//    return ostream<<query._qseqid<<'\t' <<query._sseqid<<'\t'  <<query._pident<<'\t'<<
//                    query._length<<'\t' <<query._mismatch<<'\t'<<
//                    query._gapopen<<'\t'<<query._qstart<<'\t'  <<query._qend<<'\t'<<
//                    query._sstart<<'\t' <<query._send<<'\t'    <<query._e_val<<'\t'<< query._coverage<<"\t"<<
//                    query._stitle<<'\t' <<query._species<<'\t' <<query._database_path<<'\t'<<
//                    query._frame;
}

const std::string &QuerySequence::getFrame() const {
    return _frame;
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
    _is_informative ? _sim_search_results.yes_no_inform = "Yes" : _sim_search_results.yes_no_inform = "No";
}

const std::string &QuerySequence::get_species() const {
    return _sim_search_results.species;
}

const std::string &QuerySequence::get_tax_scope() const {
    return _eggnog_results.tax_scope_readable;
}

const std::string &QuerySequence::get_contam_type() const {
    return _sim_search_results.contam_type;
}

bool QuerySequence::is_informative() const {
    return _is_informative;
}

void QuerySequence::set_contam_type(const std::string &_contam_type) {
    QuerySequence::_sim_search_results.contam_type = _contam_type;
}

void QuerySequence::set_is_database_hit(bool _is_database_hit) {
    QuerySequence::_is_database_hit = _is_database_hit;
}

// TODO move to eggnog module
void QuerySequence::set_eggnog_results(std::string seed_o, std::string seed_o_eval,
                                       std::string seed_score, std::string predicted,
                                       std::string go_terms, std::string kegg,
                                       std::string annotation_tax, std::string ogs,
                                       DatabaseHelper &database) {
    this->_eggnog_results.seed_ortholog = seed_o;
    this->_eggnog_results.seed_evalue = seed_o_eval;
    this->_eggnog_results.seed_score = seed_score;
    this->_eggnog_results.predicted_gene = predicted;
    this->_eggnog_results.ogs = ogs;
    this->_is_eggnog_hit = true;

    // Lookup/Assign Tax Scope
    if (!annotation_tax.empty()) {
        unsigned short p = (unsigned short) (annotation_tax.find("NOG"));
        if (p != std::string::npos) {
            this->_eggnog_results.tax_scope = annotation_tax.substr(0,p+3);
            this->_eggnog_results.tax_scope_readable =
                    EGGNOG_LEVELS[this->_eggnog_results.tax_scope];
        } else {
            this->_eggnog_results.tax_scope = annotation_tax;
            this->_eggnog_results.tax_scope_readable = "";
        }
    } else {
        this->_eggnog_results.tax_scope = "";
        this->_eggnog_results.tax_scope_readable = "";
    }

    // Push each GO term to a new position in vector array
    std::string temp;
    if (!go_terms.empty()) {
        std::istringstream ss(go_terms);
        while (std::getline(ss,temp,',')) {
            this->_eggnog_results.raw_go.push_back(temp);
        }
    }

    // Push each KEGG term to a new position in vector array
    if (!kegg.empty()) {
        std::istringstream keggs(kegg);
        while (std::getline(keggs,temp,',')) {
            this->_eggnog_results.raw_kegg.push_back(temp);
        }
    }

    // Find OG query was assigned to
    if (!ogs.empty()) {
        std::istringstream ss(ogs);
        std::unordered_map<std::string,std::string> og_map; // Not fully used right now
        while (std::getline(ss,temp,',')) {
            unsigned short p = (unsigned short) temp.find("@");
            og_map[temp.substr(p+1)] = temp.substr(0,p);
        }
        _eggnog_results.og_key = "";
        if (og_map.find(_eggnog_results.tax_scope) != og_map.end()) {
            _eggnog_results.og_key = og_map[_eggnog_results.tax_scope];
        }
    }

    // Lookup description, KEGG, protein domain from SQL database
    _eggnog_results.sql_kegg = "";
    _eggnog_results.description = "";
    _eggnog_results.protein_domains = "";
    if (!_eggnog_results.og_key.empty()) {
        std::vector<std::vector<std::string>>results;
        std::string sql_kegg;
        std::string sql_desc;
        std::string sql_protein;

        char *query = sqlite3_mprintf(
                "SELECT description, KEGG_freq, SMART_freq FROM og WHERE og=%Q",
                _eggnog_results.og_key.c_str());
        try {
            results = database.query(query);
            sql_desc = results[0][0];
            sql_kegg = results[0][1];
            sql_protein = results[0][2];
            if (!sql_desc.empty() && sql_desc.find("[]") != 0) _eggnog_results.description = sql_desc;
            if (!sql_kegg.empty() && sql_kegg.find("[]") != 0) {
                _eggnog_results.sql_kegg = format_eggnog(sql_kegg);
            }
            if (!sql_protein.empty() && sql_protein.find("{}") != 0){
                _eggnog_results.protein_domains = format_eggnog(sql_protein);
            }
        } catch (std::exception &e) {
            // Do not fatal error
            print_debug(e.what());
        }
    }
}


void QuerySequence::set_go_parsed(const QuerySequence::go_struct &_go_parsed, short i) {
    switch (i) {
        case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
            _eggnog_results.parsed_go = _go_parsed;
            break;
        case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
            break;
        default:
            break;
    }
}

void QuerySequence::init_sequence() {
    _seq_length = 0;
    _e_val = 0;
    _coverage = 0;

    _sim_search_results = {};
    _eggnog_results     = {};
    _interpro_results   = {};

    _frame = "";
    _sequence_p = "";
    _sequence_n = "";
    _is_family_assigned = false;
    _is_one_go = false;
    _is_one_kegg = false;
    _is_database_hit = false;
    _is_expression_kept = false;
    _is_eggnog_hit  = false;
    _is_interpro_hit = false;
}

void QuerySequence::set_lineage(const std::string &_lineage) {
    QuerySequence::_sim_search_results.lineage = _lineage;
}

void QuerySequence::set_tax_score(std::string input_lineage) {
    float tax_score = 0;
    std::string lineage = _sim_search_results.lineage;
    std::remove_if(lineage.begin(),lineage.end(), ::isspace);
    std::remove_if(input_lineage.begin(),input_lineage.end(), ::isspace);

    std::string temp;
    size_t p = 0;std::string del = ";";
    while ((p = lineage.find(";"))!=std::string::npos) {
        temp = lineage.substr(0,p);
        if (input_lineage.find(temp)!=std::string::npos) tax_score++;
        lineage.erase(0,p+del.length());
    }
    if (tax_score == 0) {
        if(_is_informative) tax_score += INFORM_ADD;
    } else {
        tax_score *= INFORM_FACTOR;
    }
    _tax_score = tax_score;
}

const std::string &QuerySequence::get_sequence() const {
    if (_sequence_n.empty()) return _sequence_p;
    return _sequence_n;
}

const QuerySequence::SimSearchResults &QuerySequence::get_sim_struct() const {
    return _sim_search_results;
}

void QuerySequence::set_sim_struct(const SimSearchResults &sim) {
    _sim_search_results = sim;
}

void QuerySequence::setIs_protein(bool is_protein) {
    QuerySequence::is_protein = is_protein;
}

bool QuerySequence::is_is_database_hit() const {
    return _is_database_hit;
}

bool QuerySequence::is_is_family_assigned() const {
    return _is_family_assigned;
}

void QuerySequence::set_is_family_assigned(bool _is_family_assigned) {
    QuerySequence::_is_family_assigned = _is_family_assigned;
}

bool QuerySequence::is_is_one_go() const {
    return _is_one_go;
}

void QuerySequence::set_is_one_go(bool _is_one_go) {
    QuerySequence::_is_one_go = _is_one_go;
}

bool QuerySequence::is_is_one_kegg() const {
    return _is_one_kegg;
}

void QuerySequence::set_is_one_kegg(bool _is_one_kegg) {
    QuerySequence::_is_one_kegg = _is_one_kegg;
}

bool QuerySequence::is_is_expression_kept() const {
    return _is_expression_kept;
}

void QuerySequence::set_is_expression_kept(bool _is_expression_kept) {
    QuerySequence::_is_expression_kept = _is_expression_kept;
};

std::string QuerySequence::print_tsv(const std::vector<const std::string*>& headers) {
    std::stringstream ss;

    // Fix, shouldn't be initialized here
    init_header();

    for (const std::string *header : headers) {
        ss << *OUTPUT_MAP[header] << "\t";
    }
    return ss.str();
}

std::string QuerySequence::print_tsv(std::vector<const std::string*>& headers, short lvl) {
    init_header();
    std::stringstream stream;
    go_struct go_terms;

    for (const std::string *header : headers) {

        if ((*header == ENTAP_EXECUTE::HEADER_EGG_GO_BIO) ||
            (*header == ENTAP_EXECUTE::HEADER_EGG_GO_CELL)||
            (*header == ENTAP_EXECUTE::HEADER_EGG_GO_MOLE)) {
            go_terms = _eggnog_results.parsed_go;
        } else go_terms = _interpro_results.parsed_go;

        if ((*header == ENTAP_EXECUTE::HEADER_EGG_GO_BIO) ||
            (*header == ENTAP_EXECUTE::HEADER_INTER_GO_BIO)) {
            if (go_terms.empty()) {
                stream <<'\t';continue;
            }
            for (std::string val : go_terms[ENTAP_EXECUTE::GO_BIOLOGICAL_FLAG]) {
                if (val.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                    stream<<val<<",";
                }
            }
            stream<<'\t';
        } else if ((*header == ENTAP_EXECUTE::HEADER_EGG_GO_CELL) ||
                   (*header == ENTAP_EXECUTE::HEADER_INTER_GO_CELL)) {
            if (go_terms.empty()) {
                stream <<'\t';continue;
            }
            for (std::string val : go_terms[ENTAP_EXECUTE::GO_CELLULAR_FLAG])  {
                if (val.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                    stream<<val<<",";
                }
            }
            stream<<'\t';
        } else if ((*header == ENTAP_EXECUTE::HEADER_EGG_GO_MOLE) ||
                   (*header == ENTAP_EXECUTE::HEADER_INTER_GO_MOLE)) {
            if (go_terms.empty()) {
                stream <<'\t';continue;
            }
            for (std::string val : go_terms[ENTAP_EXECUTE::GO_MOLECULAR_FLAG]) {
                if (val.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                    stream<<val<<",";
                }
            }
            stream<<'\t';
        } else stream << *OUTPUT_MAP[header] << "\t";
    }
    return stream.str();
}

void QuerySequence::set_fpkm(float _fpkm) {
    QuerySequence::_fpkm = _fpkm;
}


void QuerySequence::init_header() {
    OUTPUT_MAP = {
            {&ENTAP_EXECUTE::HEADER_QUERY           , &_sim_search_results.qseqid},
            {&ENTAP_EXECUTE::HEADER_SUBJECT         , &_sim_search_results.sseqid},
            {&ENTAP_EXECUTE::HEADER_PERCENT         , &_sim_search_results.pident},
            {&ENTAP_EXECUTE::HEADER_ALIGN_LEN       , &_sim_search_results.length},
            {&ENTAP_EXECUTE::HEADER_MISMATCH        , &_sim_search_results.mismatch},
            {&ENTAP_EXECUTE::HEADER_GAP_OPEN        , &_sim_search_results.gapopen},
            {&ENTAP_EXECUTE::HEADER_QUERY_E         , &_sim_search_results.qend},
            {&ENTAP_EXECUTE::HEADER_QUERY_S         , &_sim_search_results.qstart},
            {&ENTAP_EXECUTE::HEADER_SUBJ_S          , &_sim_search_results.sstart},
            {&ENTAP_EXECUTE::HEADER_SUBJ_E          , &_sim_search_results.send},
            {&ENTAP_EXECUTE::HEADER_E_VAL           , &_sim_search_results.e_val},
            {&ENTAP_EXECUTE::HEADER_COVERAGE        , &_sim_search_results.coverage},
            {&ENTAP_EXECUTE::HEADER_TITLE           , &_sim_search_results.stitle},
            {&ENTAP_EXECUTE::HEADER_SPECIES         , &_sim_search_results.species},
            {&ENTAP_EXECUTE::HEADER_DATABASE        , &_sim_search_results.database_path},
            {&ENTAP_EXECUTE::HEADER_FRAME           , &_frame},
            {&ENTAP_EXECUTE::HEADER_CONTAM          , &_sim_search_results.yes_no_contam},
            {&ENTAP_EXECUTE::HEADER_INFORM          , &_sim_search_results.yes_no_inform},
            {&ENTAP_EXECUTE::HEADER_SEED_ORTH       , &_eggnog_results.seed_ortholog},
            {&ENTAP_EXECUTE::HEADER_SEED_EVAL       , &_eggnog_results.seed_evalue},
            {&ENTAP_EXECUTE::HEADER_SEED_SCORE      , &_eggnog_results.seed_score},
            {&ENTAP_EXECUTE::HEADER_PRED_GENE       , &_eggnog_results.predicted_gene},
            {&ENTAP_EXECUTE::HEADER_TAX_SCOPE       , &_eggnog_results.tax_scope_readable},
            {&ENTAP_EXECUTE::HEADER_EGG_OGS         , &_eggnog_results.ogs},
            {&ENTAP_EXECUTE::HEADER_EGG_DESC        , &_eggnog_results.description},
            {&ENTAP_EXECUTE::HEADER_EGG_KEGG        , &_eggnog_results.sql_kegg} ,
            {&ENTAP_EXECUTE::HEADER_EGG_PROTEIN     , &_eggnog_results.protein_domains},
            {&ENTAP_EXECUTE::HEADER_INTER_EVAL      , &_interpro_results.e_value},
            {&ENTAP_EXECUTE::HEADER_INTER_INTERPRO  , &_interpro_results.interpro_desc_id},
            {&ENTAP_EXECUTE::HEADER_INTER_DATA_TERM , &_interpro_results.database_desc_id},
            {&ENTAP_EXECUTE::HEADER_INTER_DATA_TYPE , &_interpro_results.database_type},
            {&ENTAP_EXECUTE::HEADER_INTER_PATHWAY   , &_interpro_results.pathways}
    };
}


std::string QuerySequence::format_eggnog(std::string& input) {

    enum FLAGS {
        DOM    = 0x01,
        INNER  = 0x02,
        INNER2 = 0x04,
        STR    = 0x08,
        FOUND  = 0x10
    };

    unsigned char bracketFlag = 0x0;
    std::string output = "";

    for (char c : input) {
        if (c == '{') {
            bracketFlag |= DOM;
        } else if (c == '[' && !(bracketFlag & INNER)) {
            bracketFlag |= INNER;
            if (bracketFlag & DOM) output += " (";
        } else if (c == '[' && !(bracketFlag & INNER2)) {
            bracketFlag |= INNER2;
        } else if (c == '}') {
            bracketFlag &= ~DOM;
            bracketFlag &= ~FOUND;
            output = output.substr(0,output.length()-2); // Remove trailing ', '
        } else if (c == ']' && (bracketFlag & INNER2)) {
            bracketFlag &= ~INNER2;
            bracketFlag &= ~FOUND;
            output += ", ";
        } else if (c == ']' && (bracketFlag & INNER)) {
            bracketFlag &= ~INNER;
            bracketFlag &= ~FOUND;
            output = output.substr(0,output.length()-2); // Remove trailing ', '
            if (bracketFlag & DOM) output += "), ";
        } else if (c == '\"') {
            bracketFlag ^= STR;
            if (!(bracketFlag & STR)) bracketFlag |= FOUND;
            if (!(bracketFlag & INNER)) bracketFlag &= ~FOUND;
        } else {
            if (bracketFlag & FOUND) continue;
            if (bracketFlag & STR) output += c;
        }
    }
    return output;
}

void QuerySequence::set_interpro_results(std::string& eval, std::string& database_info, std::string& data,
                                         std::string& interpro_info, std::string& pathway,
                                         QuerySequence::go_struct& go_terms) {
    this->_is_interpro_hit                   = true;
    this->_interpro_results.database_desc_id = database_info;
    this->_interpro_results.database_type    = data;
    this->_interpro_results.parsed_go        = go_terms;
    this->_interpro_results.interpro_desc_id = interpro_info;
    this->_interpro_results.pathways         = pathway;
    this->_interpro_results.e_value          = eval;
}

void QuerySequence::set_is_interpro_hit(bool _is_interpro_hit) {
    QuerySequence::_is_interpro_hit = _is_interpro_hit;
}
