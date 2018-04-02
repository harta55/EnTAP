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


#include <iostream>
#include <algorithm>
#include <sstream>
#include <netinet/in.h>
#include "QuerySequence.h"
#include "EntapGlobals.h"
#include "FileSystem.h"
#include "common.h"
#include "ExceptionHandler.h"

unsigned long QuerySequence::getSeq_length() const {
    return _seq_length;
}

QuerySequence::QuerySequence() {
    init_sequence();

}

void QuerySequence::setSequence( std::string &seq) {
    QUERY_FLAG_SET(QUERY_IS_PROTEIN);
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

QuerySequence::QuerySequence(bool is_protein, std::string seq, std::string seqid){
    init_sequence();
    this->_seq_id = seqid;
    is_protein ? this->QUERY_FLAG_SET(QUERY_IS_PROTEIN) : this->QUERY_FLAG_CLEAR(QUERY_IS_PROTEIN);
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

const std::string &QuerySequence::get_species() const {
    return _sim_search_alignment_data->results->species;
}

const std::string &QuerySequence::get_contam_type() const {
    return _sim_search_alignment_data->results->contam_type;
}


void QuerySequence::set_eggnog_results(const EggnogResults &eggnogResults) {
    this->_eggnog_results = eggnogResults;
    this->QUERY_FLAG_SET(QUERY_EGGNOG_HIT);
    this->QUERY_FLAG_SET(QUERY_FAMILY_ASSIGNED);

}

void QuerySequence::init_sequence() {
    _seq_length = 0;
    _fpkm = 0;

    _sim_search_alignment_data = new SimSearchAlignmentData();
    _eggnog_results     = {};
    _interpro_results   = {};

    _frame = "";
    _sequence_p = "";
    _sequence_n = "";

    _query_flags = 0;
    QUERY_FLAG_SET(QUERY_FRAME_KEPT);
    QUERY_FLAG_SET(QUERY_EXPRESSION_KEPT);
}

const std::string &QuerySequence::get_sequence() const {
    if (_sequence_n.empty()) return _sequence_p;
    return _sequence_n;
}

std::string QuerySequence::print_tsv(const std::vector<const std::string*>& headers) {
    std::stringstream ss;

    // Fix, shouldn't be initialized here
    init_header();

    for (const std::string *header : headers) {
        ss << *OUTPUT_MAP[header] << "\t";
    }
    return ss.str();
}

/**
 * ======================================================================
 * Function std::string QuerySequence::print_tsv(std::vector<const std::string*>& headers,
 *                                                  short lvl)
 *
 * Description          - Formats data from Query Sequence to be print in ontology
 *                      -
 *
 * Notes                - None
 *
 * @param headers       - Map of headers from ontology
 * @param lvl           - Go level that would be normalized to
 *
 * @return              - String of output
 *
 * =====================================================================
 */
std::string QuerySequence::print_tsv(std::vector<const std::string*>& headers, short lvl) {
    init_header();
    std::stringstream stream;
    go_format_t go_terms;

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
        } else {
            if (OUTPUT_MAP[header] == nullptr) {
                stream << "" << "\t";
            } else {
                stream << *OUTPUT_MAP[header] << "\t";
            }
        }
    }
    // free memory in case it gets set later (unlikely)
    return stream.str();
}

void QuerySequence::set_fpkm(float _fpkm) {
    QuerySequence::_fpkm = _fpkm;
}


void QuerySequence::init_header() {

    OUTPUT_MAP = {
            {&ENTAP_EXECUTE::HEADER_QUERY           , &_seq_id},
            {&ENTAP_EXECUTE::HEADER_SUBJECT         , &_sim_search_alignment_data->results->sseqid},
            {&ENTAP_EXECUTE::HEADER_PERCENT         , &_sim_search_alignment_data->results->pident},
            {&ENTAP_EXECUTE::HEADER_ALIGN_LEN       , &_sim_search_alignment_data->results->length},
            {&ENTAP_EXECUTE::HEADER_MISMATCH        , &_sim_search_alignment_data->results->mismatch},
            {&ENTAP_EXECUTE::HEADER_GAP_OPEN        , &_sim_search_alignment_data->results->gapopen},
            {&ENTAP_EXECUTE::HEADER_QUERY_E         , &_sim_search_alignment_data->results->qend},
            {&ENTAP_EXECUTE::HEADER_QUERY_S         , &_sim_search_alignment_data->results->qstart},
            {&ENTAP_EXECUTE::HEADER_SUBJ_S          , &_sim_search_alignment_data->results->sstart},
            {&ENTAP_EXECUTE::HEADER_SUBJ_E          , &_sim_search_alignment_data->results->send},
            {&ENTAP_EXECUTE::HEADER_E_VAL           , &_sim_search_alignment_data->results->e_val},
            {&ENTAP_EXECUTE::HEADER_COVERAGE        , &_sim_search_alignment_data->results->coverage},
            {&ENTAP_EXECUTE::HEADER_TITLE           , &_sim_search_alignment_data->results->stitle},
            {&ENTAP_EXECUTE::HEADER_SPECIES         , &_sim_search_alignment_data->results->species},
            {&ENTAP_EXECUTE::HEADER_DATABASE        , &_sim_search_alignment_data->results->database_path},
            {&ENTAP_EXECUTE::HEADER_FRAME           , &_frame},
            {&ENTAP_EXECUTE::HEADER_CONTAM          , &_sim_search_alignment_data->results->yes_no_contam},
            {&ENTAP_EXECUTE::HEADER_INFORM          , &_sim_search_alignment_data->results->yes_no_inform},
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

    if (!QUERY_FLAG_GET(QUERY_BLAST_HIT)) {
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_SUBJECT] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_PERCENT] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_ALIGN_LEN] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_MISMATCH] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_GAP_OPEN] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_QUERY_E] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_QUERY_S] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_SUBJ_S] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_SUBJ_E] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_E_VAL] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_COVERAGE] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_TITLE] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_SPECIES] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_DATABASE] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_CONTAM] = nullptr;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_INFORM] = nullptr;
    }
}

void QuerySequence::set_interpro_results(std::string& eval, std::string& database_info, std::string& data,
                                         std::string& interpro_info, std::string& pathway,
                                         go_format_t& go_terms) {
    this->QUERY_FLAG_SET(QUERY_INTERPRO);
    this->_interpro_results.database_desc_id = database_info;
    this->_interpro_results.database_type    = data;
    this->_interpro_results.parsed_go        = go_terms;
    this->_interpro_results.interpro_desc_id = interpro_info;
    this->_interpro_results.pathways         = pathway;
    this->_interpro_results.e_value          = eval;
}

bool QuerySequence::is_kept() {
    return QUERY_FLAG_GET(QUERY_EXPRESSION_KEPT) &&
            QUERY_FLAG_GET(QUERY_FRAME_KEPT);
}

bool QuerySequence::QUERY_FLAG_GET(QUERY_FLAGS flag) {
    return (_query_flags & flag) != 0;
}

void QuerySequence::QUERY_FLAG_SET(QUERY_FLAGS flag) {
    _query_flags |= flag;
}

void QuerySequence::QUERY_FLAG_CLEAR(QUERY_FLAGS flag) {
    _query_flags &= ~flag;
}


QuerySequence::~QuerySequence() {
    for (auto &pair : (_sim_search_alignment_data->alignments)) {
        // first: database, second: pair of best hit to all hits (vector)
        for (QueryAlignment *alignment : pair.second.second) {
            // cycle through alignments and update the best hit
            delete(alignment);
        }
    }
    delete _sim_search_alignment_data;
}

// **********************************************************************


SimSearchAlignment::SimSearchAlignment(QuerySequence * a, uint16 b, std::string c, SimSearchResults d, std::string &lineage)
        : QueryAlignment(a,b,c){
    this->_sim_search_results = d;
    this->_best_hit = false;
    this->set_tax_score(lineage);

    this->OUTPUT_MAP = {
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
            {&ENTAP_EXECUTE::HEADER_INFORM          , &_sim_search_results.yes_no_inform}
    };
}

SimSearchResults *SimSearchAlignment::get_results() {
    return &_sim_search_results;
}

bool SimSearchAlignment::operator>(const SimSearchAlignment &alignment) {

    fp64 eval1 = this->_sim_search_results.e_val_raw;
    fp64 eval2 = alignment._sim_search_results.e_val_raw;
    // Avoid error on taking log
    if (eval1 == 0) eval1 = 1E-200;
    if (eval2 == 0) eval2 = 1E-200;
    fp64 cov1 = this->_sim_search_results.coverage_raw;
    fp64 cov2 = fabs(alignment._sim_search_results.coverage_raw);
    fp64 coverage_dif = fabs(cov1 - cov2);
    if (!this->_best_hit) {
        // For hits of the same database "better hit"
        if (fabs(log10(eval1) - log10(eval2)) < E_VAL_DIF) {
            if (coverage_dif > COV_DIF) {
                return cov1 > cov2;
            }
            if (this->_sim_search_results.contaminant && !alignment._sim_search_results.contaminant) return false;
            if (!this->_sim_search_results.contaminant && alignment._sim_search_results.contaminant) return true;
            if (this->_sim_search_results.tax_score == alignment._sim_search_results.tax_score)
                return eval1 < eval2;
            return this->_sim_search_results.tax_score > alignment._sim_search_results.tax_score;
        } else {
            return eval1 < eval2;
        }
    }else {
        // For overall best hits between databases "best hit"
        if (coverage_dif > COV_DIF) {
            return cov1 > cov2;
        }
        if (this->_sim_search_results.contaminant && !alignment._sim_search_results.contaminant) return false;
        if (!this->_sim_search_results.contaminant && alignment._sim_search_results.contaminant) return true;
        return this->_sim_search_results.tax_score > alignment._sim_search_results.tax_score;
    }
}

QueryAlignment::QueryAlignment(QuerySequence *parent, uint16 flag, std::string database) {
    _software_flag = flag;
    _parent_sequence = parent;
    _database = database;
    _frame = parent->getFrame();
}

std::string QueryAlignment::print_tsv(const std::vector<const std::string *> &headers)  {
    std::stringstream stream;

    for (const std::string *header : headers) {
        stream << *OUTPUT_MAP[header] << "\t";
    }
    return stream.str();
}

bool QuerySequence::hit_database(std::string &database, ExecuteStates state) {
    switch (state) {
        case SIMILARITY_SEARCH:
            if (_sim_search_alignment_data == nullptr) return false;
            if (_sim_search_alignment_data->alignments.empty()) return false;
            if (!database.empty()) {
                return (_sim_search_alignment_data->alignments.find(database) !=\
                    _sim_search_alignment_data->alignments.end());
            } else {
                return QUERY_FLAG_GET(QUERY_BLAST_HIT);
            }
        default:
            return false;
    }
}


QuerySequence::align_database_hits_t* QuerySequence::get_database_hits(std::string &database, ExecuteStates state) {
    switch (state) {
        case SIMILARITY_SEARCH:
            return &this->_sim_search_alignment_data->alignments[database];
        default:
            return nullptr;
    }
}

void QuerySequence::update_query_flags(ExecuteStates state) {
    switch (state) {
        case SIMILARITY_SEARCH:
            _sim_search_alignment_data-> results = _sim_search_alignment_data->best_hit->get_results();
            if (this->_sim_search_alignment_data->results->is_informative) {
                QUERY_FLAG_SET(QUERY_INFORMATIVE);
            } else {
                QUERY_FLAG_CLEAR(QUERY_INFORMATIVE);
            }
            if (this->_sim_search_alignment_data->results->contaminant) {
                QUERY_FLAG_SET(QUERY_CONTAMINANT);
            } else {
                QUERY_FLAG_CLEAR(QUERY_CONTAMINANT);
            }
            break;
        default:
            break;
    }

}

bool QuerySequence::isContaminant() {
    return this->QUERY_FLAG_GET(QUERY_CONTAMINANT);
}

/**
 * ======================================================================
 * Function void QuerySequence::set_tax_score(std::string input_lineage)
 *
 * Description          - Sets tax score based on informativeness and
 *                        lineage
 *
 * Notes                - None
 *
 * @param lineage       - Lineage input from  user
 *
 * @return              - None
 *
 * =====================================================================
 */
void SimSearchAlignment::set_tax_score(std::string &input_lineage) {
    float tax_score = 0;
    std::string lineage = this->_sim_search_results.lineage;
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
        if(this->_sim_search_results.is_informative) tax_score += INFORM_ADD;
    } else {
        if (this->_sim_search_results.is_informative) tax_score *= INFORM_FACTOR;
    }
    this->_sim_search_results.tax_score = tax_score;
}


