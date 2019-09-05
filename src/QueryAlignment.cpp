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

#include "QueryAlignment.h"


//**********************************************************************
//**********************************************************************
//                 QueryAlignment Class
//**********************************************************************
//**********************************************************************

QueryAlignment::QueryAlignment(ExecuteStates state, uint16 software, std::string &database_path, QuerySequence* parent) {
    mExecutionState = state;
    mSoftwareModule = software;
    mDatabasePath   = database_path;
    mpParentSequence= parent;
    mCompareOverallAlignment = false;

}

void QueryAlignment::set_compare_overall_alignment(bool val) {
    mCompareOverallAlignment = val;
}

std::string QueryAlignment::print_delim(std::vector<ENTAP_HEADERS> &headers, uint8 lvl, char delim)  {
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
                mpParentSequence->get_header_data(temp, header, lvl);
                stream << temp << delim;
            }
        }
    }
    return stream.str();
}

void QueryAlignment::get_all_header_data(std::string *headers) {
    for (auto &pair : ALIGN_OUTPUT_MAP) {
        headers[pair.first] = *pair.second;
    }
}

void QueryAlignment::get_header_data(ENTAP_HEADERS header, std::string &val, uint8 lvl) {
    std::vector<std::string> go_list;

    if (is_go_header(header, go_list)) {
        val = mpParentSequence->format_go_info(go_list, lvl);
    } else {
        val = *ALIGN_OUTPUT_MAP[header];
    }
}

uint16 QueryAlignment::getMSoftwareModule() const {
    return mSoftwareModule;
}

ExecuteStates QueryAlignment::getMExecutionState() const {
    return mExecutionState;
}

std::string &QueryAlignment::getMDatabasePath() {
    return mDatabasePath;
}


//**********************************************************************
//**********************************************************************
//                 SimSearchAlignment Class
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
void SimSearchAlignment::set_tax_score(std::string &input_lineage) {
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


SimSearchAlignment::SimSearchAlignment(ExecuteStates state, uint16 software, std::string &database_path, QuerySequence* parent,
                                       QuerySequence::SimSearchResults d, std::string &lineage)
    : QueryAlignment(state, software, database_path, parent){
    _sim_search_results = d;
    set_tax_score(lineage);

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
            {ENTAP_HEADER_SIM_TAXONOMIC_LINEAGE, &_sim_search_results.lineage},
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

QuerySequence::SimSearchResults* SimSearchAlignment::get_results() {
    return &_sim_search_results;
}

bool SimSearchAlignment::operator>(const QueryAlignment &alignment) {

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
    if (!this->mCompareOverallAlignment) {
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
        if (this->_sim_search_results.tax_score == alignment_cast._sim_search_results.tax_score) {
			return cov1 > cov2;
		} else {
			return this->_sim_search_results.tax_score > alignment_cast._sim_search_results.tax_score;
		}
    }
}

bool SimSearchAlignment::is_go_header(ENTAP_HEADERS header, std::vector<std::string> &go_list) {

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
//                 EggnogDmndAlignment Class
//**********************************************************************
//**********************************************************************

EggnogDmndAlignment::EggnogDmndAlignment(ExecuteStates state, uint16 software, std::string &database_path, QuerySequence* parent,
                                         QuerySequence::EggnogResults eggnogResults)
    : QueryAlignment (state, software, database_path, parent) {
    mEggnogResults = eggnogResults;
    refresh_headers();
}


QuerySequence::EggnogResults *EggnogDmndAlignment::get_results() {
    return &this->mEggnogResults;
}

bool EggnogDmndAlignment::operator>(const QueryAlignment & alignment) {
    const EggnogDmndAlignment alignment_cast = dynamic_cast<const EggnogDmndAlignment&>(alignment);

    return this->mEggnogResults.seed_eval_raw < alignment_cast.mEggnogResults.seed_eval_raw;
}

bool EggnogDmndAlignment::is_go_header(ENTAP_HEADERS header, std::vector<std::string> &go_list) {
    bool out_flag;

    switch (header) {

        case ENTAP_HEADER_ONT_EGG_GO_CELL:
            go_list = mEggnogResults.parsed_go[GO_CELLULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_ONT_EGG_GO_MOLE:
            go_list = mEggnogResults.parsed_go[GO_MOLECULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_ONT_EGG_GO_BIO:
            go_list = mEggnogResults.parsed_go[GO_BIOLOGICAL_FLAG];
            out_flag = true;
            break;

        default:
            out_flag = false;
    }
    return out_flag;
}

void EggnogDmndAlignment::refresh_headers() {

    ALIGN_OUTPUT_MAP = {
            {ENTAP_HEADER_ONT_EGG_SEED_ORTHO, &mEggnogResults.seed_ortholog},
            {ENTAP_HEADER_ONT_EGG_SEED_EVAL,  &mEggnogResults.seed_evalue},
            {ENTAP_HEADER_ONT_EGG_SEED_SCORE, &mEggnogResults.seed_score},
            {ENTAP_HEADER_ONT_EGG_PRED_GENE,  &mEggnogResults.predicted_gene},
            {ENTAP_HEADER_ONT_EGG_TAX_SCOPE_READABLE,  &mEggnogResults.tax_scope_readable},
            {ENTAP_HEADER_ONT_EGG_TAX_SCOPE_MAX, &mEggnogResults.tax_scope_lvl_max},
            {ENTAP_HEADER_ONT_EGG_MEMBER_OGS,&mEggnogResults.member_ogs},
            {ENTAP_HEADER_ONT_EGG_DESC,       &mEggnogResults.description},
            {ENTAP_HEADER_ONT_EGG_BIGG,       &mEggnogResults.bigg},
            {ENTAP_HEADER_ONT_EGG_KEGG,       &mEggnogResults.kegg},
            {ENTAP_HEADER_ONT_EGG_PROTEIN,    &mEggnogResults.protein_domains},
    };
    mpParentSequence->set_header_data();
    mpParentSequence->update_query_flags(GENE_ONTOLOGY, ONT_EGGNOG_DMND);
}

//**********************************************************************
//**********************************************************************
//                 InterproAlignment Class
//**********************************************************************
//**********************************************************************

InterproAlignment::InterproAlignment(ExecuteStates state, uint16 software, std::string &database_path, QuerySequence* parent,
                                     QuerySequence::InterProResults results)
    : QueryAlignment(state, software, database_path, parent){

    mInterproResults = results;

    ALIGN_OUTPUT_MAP = {
            {ENTAP_HEADER_ONT_INTER_EVAL, &mInterproResults.e_value},
            {ENTAP_HEADER_ONT_INTER_INTERPRO, &mInterproResults.interpro_desc_id},
            {ENTAP_HEADER_ONT_INTER_DATA_TERM,&mInterproResults.database_desc_id},
            {ENTAP_HEADER_ONT_INTER_DATA_TYPE,&mInterproResults.database_type},
            {ENTAP_HEADER_ONT_INTER_PATHWAYS, &mInterproResults.pathways}
    };
}

QuerySequence::InterProResults *InterproAlignment::get_results() {
    return &this->mInterproResults;
}

bool InterproAlignment::operator>(const QueryAlignment &alignment) {
    const InterproAlignment alignment_cast = dynamic_cast<const InterproAlignment&>(alignment);

    return this->mInterproResults.e_value_raw < alignment_cast.mInterproResults.e_value_raw;
}

bool InterproAlignment::is_go_header(ENTAP_HEADERS header, std::vector<std::string> &go_list) {
    bool out_flag;

    switch (header) {

        case ENTAP_HEADER_ONT_INTER_GO_CELL:
            go_list = mInterproResults.parsed_go[GO_CELLULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_ONT_INTER_GO_MOLE:
            go_list = mInterproResults.parsed_go[GO_MOLECULAR_FLAG];
            out_flag = true;
            break;
        case ENTAP_HEADER_ONT_INTER_GO_BIO:
            go_list = mInterproResults.parsed_go[GO_BIOLOGICAL_FLAG];
            out_flag = true;
            break;

        default:
            out_flag = false;
    }
    return out_flag;
}
