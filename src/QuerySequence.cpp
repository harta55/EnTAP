/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
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
    std::string sub = seq.substr(seq.find('n')+1);
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
    set_header_data();
}

#ifdef EGGNOG_MAPPER
void QuerySequence::set_eggnog_results(const EggnogResults &eggnogResults) {
    memcpy(&this->_eggnog_results, &eggnogResults, sizeof(eggnogResults));
    QUERY_FLAG_SET(QUERY_EGGNOG_HIT);
    QUERY_FLAG_SET(QUERY_FAMILY_ASSIGNED);
}
#endif

void QuerySequence::init_sequence() {
    _seq_length = 0;
    _fpkm = 0;

    _alignment_data = new AlignmentData(this);
    _eggnog_results = EggnogResults();
    _interpro_results = InterProResults();

    _frame = "";
    _sequence_p = "";
    _sequence_n = "";

    _query_flags = 0;
    QUERY_FLAG_SET(QUERY_FRAME_KEPT);
    QUERY_FLAG_SET(QUERY_EXPRESSION_KEPT);

    set_header_data();
}

const std::string &QuerySequence::get_sequence() const {
    if (_sequence_n.empty()) return _sequence_p;
    return _sequence_n;
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
std::string QuerySequence::print_delim(std::vector<ENTAP_HEADERS> &headers, short lvl, char delim) {
//    init_header();
    std::stringstream stream;
    go_format_t go_terms;
    std::string val;

    for (ENTAP_HEADERS &header : headers) {
        if (ENTAP_HEADER_INFO[header].print_header) {
            get_header_data(val, header, lvl);
            stream << val << delim;
        }
    }
    return stream.str();
}

void QuerySequence::get_header_data(std::string &data, ENTAP_HEADERS header, uint8 lvl) {
    QueryAlignment *align_ptr = nullptr;

    data = "";

    switch (header) {

        case ENTAP_HEADER_SIM_UNI_GO_BIO:
        case ENTAP_HEADER_SIM_UNI_GO_MOLE:
        case ENTAP_HEADER_SIM_UNI_GO_CELL:
            align_ptr = get_best_hit_alignment<SimSearchAlignment>(SIMILARITY_SEARCH, SIM_DIAMOND, "");
            break;

        case ENTAP_HEADER_ONT_EGG_GO_MOLE:
        case ENTAP_HEADER_ONT_EGG_GO_CELL:
        case ENTAP_HEADER_ONT_EGG_GO_BIO:
            align_ptr = get_best_hit_alignment<EggnogDmndAlignment>(GENE_ONTOLOGY, ONT_EGGNOG_DMND, "");
            break;

        case ENTAP_HEADER_ONT_INTER_GO_BIO:
        case ENTAP_HEADER_ONT_INTER_GO_CELL:
        case ENTAP_HEADER_ONT_INTER_GO_MOLE:
            align_ptr = get_best_hit_alignment<InterproAlignment>(GENE_ONTOLOGY, ONT_INTERPRO_SCAN, "");
            break;

        default:
            data = _header_info[header];
            break;
    }

    if (align_ptr != nullptr) {
        align_ptr->get_header_data(header, data, lvl);
    }
}

void QuerySequence::set_header_data() {
    QueryAlignment *align_ptr = nullptr;

    // General data
    _header_info[ENTAP_HEADER_QUERY] = this->_seq_id;

    // Frame Selection data
    _header_info[ENTAP_HEADER_FRAME] = this->_frame;

    // Expression Filtering data
    _header_info[ENTAP_HEADER_EXP_FPKM] = float_to_string(this->_fpkm);

    // Similarity Search data
    align_ptr = this->_alignment_data->get_best_align_ptr(SIMILARITY_SEARCH, SIM_DIAMOND, "");
    if (align_ptr != nullptr) {
        align_ptr->get_all_header_data(_header_info);
    }

    // Ontology EggNOG data
    align_ptr = this->_alignment_data->get_best_align_ptr(GENE_ONTOLOGY, ONT_EGGNOG_DMND, "");
    if (align_ptr != nullptr) {
        align_ptr->get_all_header_data(_header_info);
    }

    // Ontology InterProScan data
    align_ptr = this->_alignment_data->get_best_align_ptr(GENE_ONTOLOGY, ONT_INTERPRO_SCAN, "");
    if (align_ptr != nullptr) {
        align_ptr->get_all_header_data(_header_info);
    }
}

void QuerySequence::set_fpkm(float _fpkm) {
    QuerySequence::_fpkm = _fpkm;
    set_header_data();
}

bool QuerySequence::isContaminant() {
    return this->QUERY_FLAG_GET(QUERY_CONTAMINANT);
}

bool QuerySequence::is_kept() {
    return QUERY_FLAG_GET(QUERY_EXPRESSION_KEPT) && QUERY_FLAG_GET(QUERY_FRAME_KEPT);
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
    delete _alignment_data;
}

void QuerySequence::add_alignment(ExecuteStates state, uint16 software, EggnogResults &results, std::string& database) {
    QUERY_FLAG_SET(QUERY_EGGNOG_HIT);
    QUERY_FLAG_SET(QUERY_FAMILY_ASSIGNED);
    _alignment_data->update_best_hit(state, software, database, new EggnogDmndAlignment(results,this));
}

void QuerySequence::add_alignment(ExecuteStates state, uint16 software, SimSearchResults &results, std::string& database,std::string lineage) {
    QUERY_FLAG_SET(QUERY_BLAST_HIT);
    QueryAlignment *new_alignment = new SimSearchAlignment(results, lineage, this);
    _alignment_data->update_best_hit(state, software, database, new_alignment);
}

void QuerySequence::add_alignment(ExecuteStates state, uint16 software, QuerySequence::InterProResults &results,
                                  std::string &database) {
    QUERY_FLAG_SET(QUERY_INTERPRO);
    QueryAlignment *new_alignmet = new InterproAlignment(results, this);
    _alignment_data->update_best_hit(state, software, database, new_alignmet);
}

//**********************************************************************
//**********************************************************************
//                              AlignmentData
//**********************************************************************
//**********************************************************************

QuerySequence::AlignmentData::AlignmentData(QuerySequence *sequence){
    querySequence = sequence;
}

QuerySequence::AlignmentData::~AlignmentData() {
    // remove sim search alignments
    for (ALIGNMENT_DATA_T &software_data : sim_search_data) {
        // Cycle through each software data struct
        for (auto &pair : software_data) {
            // for each database, delete vector
            for (auto &alignment : pair.second) {
                delete alignment;
            }
        }
    }

    // remove ontology alignments
    for (ALIGNMENT_DATA_T &software_data : ontology_data) {
        // Cycle through each software data struct
        for (auto &pair : software_data) {
            // for each database, delete vector
            for (auto &alignment : pair.second) {
                delete alignment;
            }
        }
    }
}


void QuerySequence::AlignmentData::update_best_hit(ExecuteStates state, uint16 software, std::string &database, QueryAlignment* new_alignment) {
    ALIGNMENT_DATA_T* alignment_arr = get_software_ptr(state, software);


    // Did we hit against this database yet
    if (!hit_database(state, software, database)) {
        // No, create new vector for that database and add as best hit for database
        align_database_hits_t vect = {new_alignment};
        alignment_arr->emplace(database, vect);
    } else {
        // Yes, add alignment to list then update
        alignment_arr->at(database).push_back(new_alignment);
    }

    // Always will have hit this database, get alignments
    align_database_hits_t *database_data = &alignment_arr->at(database);

    // Sort alignments for that database (0 index is best hit)
    std::sort(database_data->begin(), database_data->end(), sort_descending_database());

    // See if this alignment is better than the overall alignment
    new_alignment->set_compare_overall_alignment(true);
    QueryAlignment* best_alignment = get_best_align_ptr(state, software, "");

    if (best_alignment == nullptr) {
        // Replace if so
        set_best_alignment(state, software, new_alignment);
    } else {
        best_alignment->set_compare_overall_alignment(true);
        if (new_alignment > best_alignment) {
            set_best_alignment(state, software, new_alignment);
        }
    }

    // Update any overall flags that may have changed with best hit changes
    querySequence->update_query_flags(state, software);
}

bool QuerySequence::AlignmentData::hit_database(ExecuteStates state, uint16 software, std::string &database) {
    // If we want overall best alignment
    if (database.empty()) {
        return get_best_align_ptr(state, software, database) != nullptr;
    } else {
        // If we hit a specific database
        return get_database_ptr(state, software, database) != nullptr;
    }
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
                case ONT_EGGNOG_DMND: {
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
    set_header_data();
}

QuerySequence::align_database_hits_t *
QuerySequence::get_database_hits(std::string &database, ExecuteStates state, uint16 software) {
    return this->_alignment_data->get_database_ptr(state, software, database);
}

std::string QuerySequence::format_go_info(std::vector<std::string> &go_list, uint8 lvl) {
    std::stringstream out;

    for (std::string &val : go_list)  {
        if (val.find("(L=" + std::to_string(lvl))!=std::string::npos || lvl == 0) {
            out<<val<<",";
        }
    }
    return out.str();
}

bool QuerySequence::hit_database(ExecuteStates state, uint16 software, std::string database) {
    return _alignment_data->hit_database(state, software, database);
}

QuerySequence::align_database_hits_t* QuerySequence::AlignmentData::get_database_ptr(ExecuteStates state, uint16 software, std::string& database) {
    if (database.empty()) return nullptr;

    switch (state) {
        case SIMILARITY_SEARCH:
            if (this->sim_search_data[software].find(database) != (this->sim_search_data[software].end())) {
                return &this->sim_search_data[software].at(database);
            } else {
                return nullptr;
            }
        case GENE_ONTOLOGY:
            if (this->ontology_data[software].find(database) != (this->ontology_data[software].end())) {
                return &this->ontology_data[software].at(database);
            } else {
                return nullptr;
            }
        default:
            return nullptr;
    }
}

QuerySequence::ALIGNMENT_DATA_T* QuerySequence::AlignmentData::get_software_ptr(ExecuteStates state, uint16 software) {
    switch (state) {
        case SIMILARITY_SEARCH:
            return &this->sim_search_data[software];
        case GENE_ONTOLOGY:
            return &this->ontology_data[software];
        default:
            return nullptr;
    }
}

QuerySequence::QueryAlignment* QuerySequence::AlignmentData::get_best_align_ptr(ExecuteStates state, uint16 software, std::string database) {
    if (database.empty()) {
        return overall_alignment[state][software];
    } else {
        return get_database_ptr(state, software, database)->at(0);
    }
}

void
QuerySequence::AlignmentData::set_best_alignment(ExecuteStates state, uint16 software, QueryAlignment *alignment) {
    overall_alignment[state][software] = alignment;
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

std::string QuerySequence::QueryAlignment::print_delim(std::vector<ENTAP_HEADERS> &headers, uint8 lvl, char delim)  {
    std::stringstream stream;
    std::string temp;

    for (ENTAP_HEADERS header : headers) {
        if (ENTAP_HEADER_INFO[header].print_header) {
            if (ALIGN_OUTPUT_MAP.find(header) != ALIGN_OUTPUT_MAP.end()) {
                // Header applies to this alignment
                get_header_data(header, temp, lvl);
                stream << temp << delim;

            } else {
                // Header does NOT apply to this alignment, get info from parent
                _parent->get_header_data(temp, header, lvl);
                stream << temp << delim;
            }
        }
    }
    return stream.str();
}

void QuerySequence::QueryAlignment::get_all_header_data(std::string *headers) {
    for (auto &pair : ALIGN_OUTPUT_MAP) {
        headers[pair.first] = *pair.second;
    }
}

void QuerySequence::QueryAlignment::get_header_data(ENTAP_HEADERS header, std::string &val, uint8 lvl) {
    std::vector<std::string> go_list;

    if (is_go_header(header, go_list)) {
        val = _parent->format_go_info(go_list, lvl);
    } else {
        val = *ALIGN_OUTPUT_MAP[header];
    }
}


//**********************************************************************
//**********************************************************************
//                 SimSearchAlignment Nested Class
//**********************************************************************
//**********************************************************************


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
    while ((p = lineage.find(';'))!=std::string::npos) {
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


QuerySequence::SimSearchAlignment::SimSearchAlignment(SimSearchResults d, std::string &lineage, QuerySequence* parent) {
    _sim_search_results = d;
    set_tax_score(lineage);
    _parent = parent;

    ALIGN_OUTPUT_MAP = {
            {ENTAP_HEADER_QUERY               , &_sim_search_results.qseqid},
            {ENTAP_HEADER_SIM_SUBJECT         , &_sim_search_results.sseqid},
            {ENTAP_HEADER_SIM_PERCENT         , &_sim_search_results.pident},
            {ENTAP_HEADER_SIM_ALIGN_LEN       , &_sim_search_results.length},
            {ENTAP_HEADER_SIM_MISMATCH        , &_sim_search_results.mismatch},
            {ENTAP_HEADER_SIM_GAP_OPEN        , &_sim_search_results.gapopen},
            {ENTAP_HEADER_SIM_QUERY_E         , &_sim_search_results.qend},
            {ENTAP_HEADER_SIM_QUERY_S         , &_sim_search_results.qstart},
            {ENTAP_HEADER_SIM_SUBJ_S          , &_sim_search_results.sstart},
            {ENTAP_HEADER_SIM_SUBJ_E          , &_sim_search_results.send},
            {ENTAP_HEADER_SIM_E_VAL           , &_sim_search_results.e_val},
            {ENTAP_HEADER_SIM_COVERAGE        , &_sim_search_results.coverage},
            {ENTAP_HEADER_SIM_TITLE           , &_sim_search_results.stitle},
            {ENTAP_HEADER_SIM_SPECIES         , &_sim_search_results.species},
            {ENTAP_HEADER_SIM_DATABASE        , &_sim_search_results.database_path},
            {ENTAP_HEADER_SIM_CONTAM          , &_sim_search_results.yes_no_contam},
            {ENTAP_HEADER_SIM_INFORM          , &_sim_search_results.yes_no_inform},
            {ENTAP_HEADER_SIM_UNI_DATA_XREF   , &_sim_search_results.uniprot_info.database_x_refs},
            {ENTAP_HEADER_SIM_UNI_COMMENTS    , &_sim_search_results.uniprot_info.comments},
//            {ENTAP_HEADER_SIM_UNI_GO_CELL     , &_sim_search_results.uniprot_info.go_terms.at(GO_CELLULAR_FLAG)},
//            {ENTAP_HEADER_SIM_UNI_GO_BIO      , &_sim_search_results.uniprot_info.go_terms.at(GO_BIOLOGICAL_FLAG)},
//            {ENTAP_HEADER_SIM_UNI_GO_MOLE     , &_sim_search_results.uniprot_info.go_terms.at(GO_MOLECULAR_FLAG)}
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

bool QuerySequence::SimSearchAlignment::is_go_header(ENTAP_HEADERS header, std::vector<std::string> &go_list) {

    bool out_flag;

    switch (header) {

        case ENTAP_HEADER_SIM_UNI_GO_CELL:
            go_list = _sim_search_results.uniprot_info.go_terms[GO_CELLULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_SIM_UNI_GO_MOLE:
            go_list = _sim_search_results.uniprot_info.go_terms[GO_MOLECULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_SIM_UNI_GO_BIO:
            go_list = _sim_search_results.uniprot_info.go_terms[GO_BIOLOGICAL_FLAG];
            out_flag = true;
            break;

        default:
            out_flag = false;
    }
    return out_flag;
}




//**********************************************************************
//**********************************************************************
//                 EggnogDmndAlignment Struct
//**********************************************************************
//**********************************************************************

QuerySequence::EggnogDmndAlignment::EggnogDmndAlignment(QuerySequence::EggnogResults eggnogResults,
                                                        QuerySequence *parent) {
    _parent = parent;
    _eggnog_results = eggnogResults;

    ALIGN_OUTPUT_MAP = {
            {ENTAP_HEADER_ONT_EGG_SEED_ORTHO, &_eggnog_results.seed_ortholog},
            {ENTAP_HEADER_ONT_EGG_SEED_EVAL,  &_eggnog_results.seed_evalue},
            {ENTAP_HEADER_ONT_EGG_SEED_SCORE, &_eggnog_results.seed_score},
            {ENTAP_HEADER_ONT_EGG_PRED_GENE,  &_eggnog_results.predicted_gene},
            {ENTAP_HEADER_ONT_EGG_TAX_SCOPE,  &_eggnog_results.tax_scope},
            {ENTAP_HEADER_ONT_EGG_MEMBER_OGS,&_eggnog_results.member_ogs},
            {ENTAP_HEADER_ONT_EGG_DESC,       &_eggnog_results.description},
            {ENTAP_HEADER_ONT_EGG_KEGG,       &_eggnog_results.kegg},
            {ENTAP_HEADER_ONT_EGG_PROTEIN,    &_eggnog_results.protein_domains},
//            {ENTAP_HEADER_ONT_EGG_GO_BIO,     &_eggnog_results.parsed_go.at(GO_BIOLOGICAL_FLAG)},
//            {ENTAP_HEADER_ONT_EGG_GO_CELL,    &_eggnog_results.parsed_go.at(GO_CELLULAR_FLAG)},
//            {ENTAP_HEADER_ONT_EGG_GO_MOLE,    &_eggnog_results.parsed_go.at(GO_MOLECULAR_FLAG)}
    };
}


QuerySequence::EggnogResults *QuerySequence::EggnogDmndAlignment::get_results() {
    return &this->_eggnog_results;
}

bool QuerySequence::EggnogDmndAlignment::operator>(const QueryAlignment & alignment) {
    const EggnogDmndAlignment alignment_cast = dynamic_cast<const EggnogDmndAlignment&>(alignment);

    return this->_eggnog_results.seed_eval_raw < alignment_cast._eggnog_results.seed_eval_raw;
}

bool QuerySequence::EggnogDmndAlignment::is_go_header(ENTAP_HEADERS header, std::vector<std::string> &go_list) {
    bool out_flag;

    switch (header) {

        case ENTAP_HEADER_ONT_EGG_GO_CELL:
            go_list = _eggnog_results.parsed_go[GO_CELLULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_ONT_EGG_GO_MOLE:
            go_list = _eggnog_results.parsed_go[GO_MOLECULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_ONT_EGG_GO_BIO:
            go_list = _eggnog_results.parsed_go[GO_BIOLOGICAL_FLAG];
            out_flag = true;
            break;

        default:
            out_flag = false;
    }
    return out_flag;
}

//**********************************************************************
//**********************************************************************
//                 InterproAlignment Struct
//**********************************************************************
//**********************************************************************

QuerySequence::InterproAlignment::InterproAlignment(QuerySequence::InterProResults results, QuerySequence *parent) {

    _interpro_results = results;
    _parent = parent;

    ALIGN_OUTPUT_MAP = {
            {ENTAP_HEADER_ONT_INTER_EVAL, &_interpro_results.e_value},
            {ENTAP_HEADER_ONT_INTER_INTERPRO, &_interpro_results.interpro_desc_id},
            {ENTAP_HEADER_ONT_INTER_DATA_TERM,&_interpro_results.database_desc_id},
            {ENTAP_HEADER_ONT_INTER_DATA_TYPE,&_interpro_results.database_type},
            {ENTAP_HEADER_ONT_INTER_PATHWAYS, &_interpro_results.pathways}
    };
}

QuerySequence::InterProResults *QuerySequence::InterproAlignment::get_results() {
    return &this->_interpro_results;
}

bool QuerySequence::InterproAlignment::operator>(const QuerySequence::QueryAlignment &alignment) {
    const InterproAlignment alignment_cast = dynamic_cast<const InterproAlignment&>(alignment);

    return this->_interpro_results.e_value_raw < alignment_cast._interpro_results.e_value_raw;
}

bool QuerySequence::InterproAlignment::is_go_header(ENTAP_HEADERS header, std::vector<std::string> &go_list) {
    bool out_flag;

    switch (header) {

        case ENTAP_HEADER_ONT_INTER_GO_CELL:
            go_list = _interpro_results.parsed_go[GO_CELLULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_ONT_INTER_GO_MOLE:
            go_list = _interpro_results.parsed_go[GO_MOLECULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_ONT_INTER_GO_BIO:
            go_list = _interpro_results.parsed_go[GO_BIOLOGICAL_FLAG];
            out_flag = true;
            break;

        default:
            out_flag = false;
    }
    return out_flag;
}