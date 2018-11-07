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

const std::string &QuerySequence::getFrame() const {
    return _frame;
}

void QuerySequence::setFrame(const std::string &frame) {
    QuerySequence::_frame = frame;
}

void QuerySequence::set_eggnog_results(const EggnogResults &eggnogResults) {
    memcpy(&this->_eggnog_results, &eggnogResults, sizeof(eggnogResults));
    this->QUERY_FLAG_SET(QUERY_EGGNOG_HIT);
    this->QUERY_FLAG_SET(QUERY_FAMILY_ASSIGNED);
}

void QuerySequence::init_sequence() {
    _seq_length = 0;
    _fpkm = 0;

    _total_alignment_data = AlignmentData();

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

        // Determine which go term data we want to use
        if ((*header == ENTAP_EXECUTE::HEADER_EGG_GO_BIO) ||
            (*header == ENTAP_EXECUTE::HEADER_EGG_GO_CELL)||
            (*header == ENTAP_EXECUTE::HEADER_EGG_GO_MOLE)) {
            go_terms = _eggnog_results.parsed_go;
        } else if ((*header == ENTAP_EXECUTE::HEADER_INTER_GO_BIO) ||
                   (*header == ENTAP_EXECUTE::HEADER_INTER_GO_CELL)||
                   (*header == ENTAP_EXECUTE::HEADER_INTER_GO_MOLE)) {
            go_terms = _interpro_results.parsed_go;
        } else if ((*header == ENTAP_EXECUTE::HEADER_UNI_GO_MOLE) ||
                   (*header == ENTAP_EXECUTE::HEADER_UNI_GO_CELL) ||
                   (*header == ENTAP_EXECUTE::HEADER_UNI_GO_BIO)){
            go_terms =  get_best_hit_alignment<SimSearchAlignment>(SIMILARITY_SEARCH,
                    ENTAP_EXECUTE::SIM_SEARCH_FLAG_DIAMOND, "")->get_results()->uniprot_info.go_terms;
        }   // cleanup!!


        if ((*header == ENTAP_EXECUTE::HEADER_EGG_GO_BIO)   ||
            (*header == ENTAP_EXECUTE::HEADER_INTER_GO_BIO) ||
            (*header == ENTAP_EXECUTE::HEADER_UNI_GO_BIO)) {
            if (go_terms.empty()) {
                stream <<'\t';continue;
            }
            for (std::string &val : go_terms[ENTAP_EXECUTE::GO_BIOLOGICAL_FLAG]) {
                if (val.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                    stream<<val<<",";
                }
            }
            stream<<'\t';
        } else if ((*header == ENTAP_EXECUTE::HEADER_EGG_GO_CELL)   ||
                   (*header == ENTAP_EXECUTE::HEADER_INTER_GO_CELL) ||
                   (*header == ENTAP_EXECUTE::HEADER_UNI_GO_CELL)) {
            if (go_terms.empty()) {
                stream <<'\t';continue;
            }
            for (std::string &val : go_terms[ENTAP_EXECUTE::GO_CELLULAR_FLAG])  {
                if (val.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
                    stream<<val<<",";
                }
            }
            stream<<'\t';
        } else if ((*header == ENTAP_EXECUTE::HEADER_EGG_GO_MOLE)   ||
                   (*header == ENTAP_EXECUTE::HEADER_INTER_GO_MOLE) ||
                   (*header == ENTAP_EXECUTE::HEADER_UNI_GO_MOLE)) {
            if (go_terms.empty()) {
                stream <<'\t';continue;
            }
            for (std::string &val : go_terms[ENTAP_EXECUTE::GO_MOLECULAR_FLAG]) {
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
    return stream.str();
}

void QuerySequence::set_fpkm(float _fpkm) {
    QuerySequence::_fpkm = _fpkm;
}


void QuerySequence::init_header() {
    OUTPUT_MAP = {
            {&ENTAP_EXECUTE::HEADER_QUERY           , &_seq_id},
            {&ENTAP_EXECUTE::HEADER_FRAME           , &_frame},
            {&ENTAP_EXECUTE::HEADER_SEED_ORTH       , &_eggnog_results.seed_ortholog},
            {&ENTAP_EXECUTE::HEADER_SEED_EVAL       , &_eggnog_results.seed_evalue},
            {&ENTAP_EXECUTE::HEADER_SEED_SCORE      , &_eggnog_results.seed_score},
            {&ENTAP_EXECUTE::HEADER_PRED_GENE       , &_eggnog_results.predicted_gene},
            {&ENTAP_EXECUTE::HEADER_TAX_SCOPE       , &_eggnog_results.tax_scope_readable},
            {&ENTAP_EXECUTE::HEADER_EGG_OGS         , &_eggnog_results.ogs},
            {&ENTAP_EXECUTE::HEADER_EGG_DESC        , &_eggnog_results.description},
            {&ENTAP_EXECUTE::HEADER_EGG_KEGG        , &_eggnog_results.kegg} ,
            {&ENTAP_EXECUTE::HEADER_EGG_PROTEIN     , &_eggnog_results.protein_domains},
            {&ENTAP_EXECUTE::HEADER_INTER_EVAL      , &_interpro_results.e_value},
            {&ENTAP_EXECUTE::HEADER_INTER_INTERPRO  , &_interpro_results.interpro_desc_id},
            {&ENTAP_EXECUTE::HEADER_INTER_DATA_TERM , &_interpro_results.database_desc_id},
            {&ENTAP_EXECUTE::HEADER_INTER_DATA_TYPE , &_interpro_results.database_type},
            {&ENTAP_EXECUTE::HEADER_INTER_PATHWAY   , &_interpro_results.pathways}
    };

    SimSearchAlignment *alignment = get_best_hit_alignment<SimSearchAlignment>(SIMILARITY_SEARCH,
                                                                               ENTAP_EXECUTE::SIM_SEARCH_FLAG_DIAMOND, "");
    if (alignment != nullptr) {
        SimSearchResults *sim_search_results = alignment->get_results();
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_SUBJECT  ] = &sim_search_results->sseqid;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_PERCENT  ] = &sim_search_results->pident;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_ALIGN_LEN] = &sim_search_results->length;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_MISMATCH ] = &sim_search_results->mismatch;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_GAP_OPEN ] = &sim_search_results->gapopen;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_QUERY_E  ] = &sim_search_results->qend;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_QUERY_S  ] = &sim_search_results->qstart;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_SUBJ_S   ] = &sim_search_results->sstart;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_SUBJ_E   ] = &sim_search_results->send;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_E_VAL    ] = &sim_search_results->e_val;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_COVERAGE ] = &sim_search_results->coverage;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_TITLE    ] = &sim_search_results->stitle;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_SPECIES  ] = &sim_search_results->species;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_DATABASE ] = &sim_search_results->database_path;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_CONTAM   ] = &sim_search_results->yes_no_contam;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_INFORM   ] = &sim_search_results->yes_no_inform;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_UNI_DATA_XREF] = &sim_search_results->uniprot_info.database_x_refs;
        OUTPUT_MAP[&ENTAP_EXECUTE::HEADER_UNI_COMMENTS]  = &sim_search_results->uniprot_info.comments;
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

void QuerySequence::QUERY_FLAG_CHANGE(QuerySequence::QUERY_FLAGS flag, bool val) {
    if (val) {
        QUERY_FLAG_SET(flag);
    } else {
        QUERY_FLAG_CLEAR(flag);
    }
}


QuerySequence::~QuerySequence() {
    // Clear alignment data
    // Cycle through vector
    for (ALIGNMENT_DATA_T &alignment_data : _total_alignment_data._alignment_data) {
        // Cycle through databases
        for (auto &pair : alignment_data) {
            // Through alignments
            for (QueryAlignment *alignment : pair.second.second) {
                delete alignment;
            }
        }
    }
}

bool QuerySequence::hit_database(ExecuteStates state, uint16 software, std::string database) {
    // If we want overall best alignment
    if (database.empty()) {
        return _total_alignment_data.index_best_align(state,software) != nullptr;
    } else {
        // If we hit a specific database
        ALIGNMENT_DATA_T *data = _total_alignment_data.index_data(state,software);
        return data->find(database) != data->end();
    }
}

void QuerySequence::update_best_hit(ExecuteStates state, uint16 software, std::string &database, QueryAlignment* new_alignment) {
    // Did we not hit against the database yet?
    if (!hit_database(state, software, database)) {
        // No, create new vector for that database and add as best hit for database
        std::vector<QueryAlignment*> vect = {new_alignment};
        _total_alignment_data.index_data(state,software)->emplace(database, std::make_pair(new_alignment, vect));
    } else {
        // Yes, add alignment to list then update
        _total_alignment_data.index_data(state,software)->at(database).second.push_back(new_alignment);
    }

    // Always will have hit this database
    align_database_hits_t *database_data =
            &_total_alignment_data.index_data(state,software)->at(database);

    // Cycle through all hits of this database
    for (auto *alignment : database_data->second) {
        // See if we should update this databases alignment
        alignment->set_compare_overall_alignment(false);
        if (database_data->first == nullptr || *alignment > *database_data->first) {
            // If it is, set that to current alignment
            database_data->first = alignment;
        }

        // See if this alignment is better than the overall alignment
        alignment->set_compare_overall_alignment(true);
        if (_total_alignment_data.index_best_align(state,software) == nullptr ||
                *alignment > *_total_alignment_data.index_best_align(state,software)) {
            _total_alignment_data.set_best_align(state,software,alignment);
        }
    }
    // Update any overall flags that may have changed with best hit changes
    update_query_flags(state, software);
}

void QuerySequence::update_query_flags(ExecuteStates state, uint16 software) {
    switch (state) {
        case SIMILARITY_SEARCH: {
            SimSearchAlignment *best_align = get_best_hit_alignment<SimSearchAlignment>(state, software, "");
            SimSearchResults *results = best_align->get_results();
            QUERY_FLAG_CHANGE(QUERY_INFORMATIVE, results->is_informative);
            QUERY_FLAG_CHANGE(QUERY_CONTAMINANT, results->contaminant);
            break;
        }

        case GENE_ONTOLOGY: {
            EggnogResults *results=nullptr;
            switch (software) {
                case ENTAP_EXECUTE::EGGNOG_DMND_INT_FLAG: {
                    EggnogDmndAlignment *best_align = get_best_hit_alignment<EggnogDmndAlignment>(state, software,"");
                    results = best_align->get_results();
                    break;
                }
                default:
                    return;
            }
            QUERY_FLAG_CHANGE(QUERY_FAMILY_ONE_GO, !results->parsed_go.empty());
            QUERY_FLAG_CHANGE(QUERY_FAMILY_ONE_KEGG, !results->kegg.empty());
        }
        default:
            break;
    }

}

bool QuerySequence::isContaminant() {
    return this->QUERY_FLAG_GET(QUERY_CONTAMINANT);
}

QuerySequence::align_database_hits_t *QuerySequence::get_database_hits(std::string &database, ExecuteStates state,uint16 software) {
    return &_total_alignment_data.index_data(state,software)->at(database);
}

std::string QuerySequence::get_header_data(const std::string *header) {
    init_header();
    return *OUTPUT_MAP[header];
}


//**********************************************************************
//**********************************************************************
//                 QueryAlignment Nested Class
//**********************************************************************
//**********************************************************************


QuerySequence::QueryAlignment::QueryAlignment() {
    _compare_overall_alignment = false;
}

void QuerySequence::QueryAlignment::set_compare_overall_alignment(bool val) {
    _compare_overall_alignment = val;
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
void QuerySequence::SimSearchAlignment::set_tax_score(std::string &input_lineage) {
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

QuerySequence::SimSearchAlignment::~SimSearchAlignment() {

}



QuerySequence::SimSearchAlignment::SimSearchAlignment(SimSearchResults d, std::string &lineage, QuerySequence* parent) {
    _sim_search_results = d;
    set_tax_score(lineage);
    _parent = parent;

    ALIGN_OUTPUT_MAP = {
            {&ENTAP_EXECUTE::HEADER_QUERY           , _sim_search_results.qseqid},
            {&ENTAP_EXECUTE::HEADER_SUBJECT         , _sim_search_results.sseqid},
            {&ENTAP_EXECUTE::HEADER_PERCENT         , _sim_search_results.pident},
            {&ENTAP_EXECUTE::HEADER_ALIGN_LEN       , _sim_search_results.length},
            {&ENTAP_EXECUTE::HEADER_MISMATCH        , _sim_search_results.mismatch},
            {&ENTAP_EXECUTE::HEADER_GAP_OPEN        , _sim_search_results.gapopen},
            {&ENTAP_EXECUTE::HEADER_QUERY_E         , _sim_search_results.qend},
            {&ENTAP_EXECUTE::HEADER_QUERY_S         , _sim_search_results.qstart},
            {&ENTAP_EXECUTE::HEADER_SUBJ_S          , _sim_search_results.sstart},
            {&ENTAP_EXECUTE::HEADER_SUBJ_E          , _sim_search_results.send},
            {&ENTAP_EXECUTE::HEADER_E_VAL           , _sim_search_results.e_val},
            {&ENTAP_EXECUTE::HEADER_COVERAGE        , _sim_search_results.coverage},
            {&ENTAP_EXECUTE::HEADER_TITLE           , _sim_search_results.stitle},
            {&ENTAP_EXECUTE::HEADER_SPECIES         , _sim_search_results.species},
            {&ENTAP_EXECUTE::HEADER_DATABASE        , _sim_search_results.database_path},
            {&ENTAP_EXECUTE::HEADER_CONTAM          , _sim_search_results.yes_no_contam},
            {&ENTAP_EXECUTE::HEADER_INFORM          , _sim_search_results.yes_no_inform},
            {&ENTAP_EXECUTE::HEADER_UNI_DATA_XREF   , _sim_search_results.uniprot_info.database_x_refs},
            {&ENTAP_EXECUTE::HEADER_UNI_COMMENTS    , _sim_search_results.uniprot_info.comments}
    };
}

QuerySequence::SimSearchResults* QuerySequence::SimSearchAlignment::get_results() {
    return &_sim_search_results;
}

bool QuerySequence::SimSearchAlignment::operator>(const QueryAlignment &alignment) {

    // Don't need to check typeid
    const SimSearchAlignment alignment_cast = dynamic_cast<const SimSearchAlignment&>(alignment);

    fp64 eval1 = this->_sim_search_results.e_val_raw;
    fp64 eval2 = alignment_cast._sim_search_results.e_val_raw;
    // Avoid error on taking log
    if (eval1 == 0) eval1 = 1E-200;
    if (eval2 == 0) eval2 = 1E-200;
    fp64 cov1 = this->_sim_search_results.coverage_raw;
    fp64 cov2 = alignment_cast._sim_search_results.coverage_raw;
    fp64 coverage_dif = fabs(cov1 - cov2);
    if (!this->_compare_overall_alignment) {
        // For hits of the same database "better hit"
        if (fabs(log10(eval1) - log10(eval2)) < E_VAL_DIF) {
            if (coverage_dif > COV_DIF) {
                return cov1 > cov2;
            }
            if (this->_sim_search_results.contaminant && !alignment_cast._sim_search_results.contaminant) return false;
            if (!this->_sim_search_results.contaminant && alignment_cast._sim_search_results.contaminant) return true;
            if (this->_sim_search_results.tax_score == alignment_cast._sim_search_results.tax_score)
                return eval1 < eval2;
            return this->_sim_search_results.tax_score > alignment_cast._sim_search_results.tax_score;
        } else {
            return eval1 < eval2;
        }
    }else {
        // For overall best hits between databases "best hit"
        if (coverage_dif > COV_DIF) {
            return cov1 > cov2;
        }
        if (this->_sim_search_results.contaminant && !alignment_cast._sim_search_results.contaminant) return false;
        if (!this->_sim_search_results.contaminant && alignment_cast._sim_search_results.contaminant) return true;
        return this->_sim_search_results.tax_score > alignment_cast._sim_search_results.tax_score;
    }
}

std::string QuerySequence::QueryAlignment::print_tsv(const std::vector<const std::string *> &headers)  {
    std::stringstream stream;

    for (const std::string *header : headers) {
        if (ALIGN_OUTPUT_MAP.find(header) != ALIGN_OUTPUT_MAP.end()) {
            stream << ALIGN_OUTPUT_MAP[header] << "\t";
        } else {
            stream << _parent->get_header_data(header) << "\t";
        }
    }
    return stream.str();
}

//**********************************************************************
//**********************************************************************
//                 AlignmentData Struct
//**********************************************************************
//**********************************************************************

QuerySequence::AlignmentData::AlignmentData() {
    _alignment_data = std::vector<ALIGNMENT_DATA_T>(ENTAP_EXECUTE::SOFTWARE_MAX * EXECUTION_MAX, ALIGNMENT_DATA_T());
    _alignment_best = std::vector<QueryAlignment*>(ENTAP_EXECUTE::SOFTWARE_MAX * EXECUTION_MAX, nullptr);
}

QuerySequence::AlignmentData::~AlignmentData() {

}

QuerySequence::ALIGNMENT_DATA_T* QuerySequence::AlignmentData::index_data(ExecuteStates state, uint16 software) {
    return &_alignment_data.at(state * ENTAP_EXECUTE::SOFTWARE_MAX + software);
}

QuerySequence::QueryAlignment *QuerySequence::AlignmentData::index_best_align(ExecuteStates state, uint16 software) {
    return _alignment_best.at(state * ENTAP_EXECUTE::SOFTWARE_MAX + software);
}

void
QuerySequence::AlignmentData::set_best_align(ExecuteStates state, uint16 software, QuerySequence::QueryAlignment *alignment) {
    _alignment_best.at(state * ENTAP_EXECUTE::SOFTWARE_MAX + software) = alignment;
}

//**********************************************************************
//**********************************************************************
//                 EggnogDmndAlignment Struct
//**********************************************************************
//**********************************************************************

QuerySequence::EggnogDmndAlignment::EggnogDmndAlignment(QuerySequence::EggnogResults eggnogResults,
                                                        QuerySequence *parent) {

}

QuerySequence::EggnogDmndAlignment::~EggnogDmndAlignment() { }


QuerySequence::EggnogResults *QuerySequence::EggnogDmndAlignment::get_results() {
    return &this->_eggnog_results;
}

bool QuerySequence::EggnogDmndAlignment::operator>(const QuerySequence::QueryAlignment & alignment) {
    const EggnogDmndAlignment alignment_cast = dynamic_cast<const EggnogDmndAlignment&>(alignment);

    return this->_eggnog_results.seed_eval_raw < alignment_cast._eggnog_results.seed_eval_raw;
}
