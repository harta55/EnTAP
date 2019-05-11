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
#include "QueryAlignment.h"

unsigned long QuerySequence::getSeq_length() const {
    return _seq_length;
}

QuerySequence::QuerySequence() {
    init_sequence();
}

const std::string &QuerySequence::get_sequence_p() const {
    return _sequence_p;
}

void QuerySequence::set_sequence_p(std::string &seq) {
    QUERY_FLAG_SET(QUERY_IS_PROTEIN);
    _seq_length = calc_seq_length(seq,true);
    if (!seq.empty() && seq[seq.length()-1] == '\n') {
        seq.pop_back();
    }
    _sequence_p = seq;}

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
    set_header_data();
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
        if (*new_alignment > *best_alignment) {
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
            switch (software) {
                case ONT_EGGNOG_DMND: {
                    EggnogDmndAlignment *best_align = get_best_hit_alignment<EggnogDmndAlignment>(state, software,"");
                    EggnogResults *results = best_align->get_results();

                    // in case results were 'refreshed'
                    if (best_align != nullptr) {
                        QUERY_FLAG_CHANGE(QUERY_FAMILY_ONE_GO, !results->parsed_go.empty());
                        QUERY_FLAG_CHANGE(QUERY_FAMILY_ONE_KEGG, !results->kegg.empty());

                        if (!results->parsed_go.empty())
                            QUERY_FLAG_SET(QUERY_ONE_GO);
                        if (!results->kegg.empty())
                            QUERY_FLAG_SET(QUERY_ONE_KEGG);
                    }
                    break;
                }
                case ONT_INTERPRO_SCAN: {
                    InterproAlignment *best_align = get_best_hit_alignment<InterproAlignment>(state, software, "");
                    InterProResults *results = best_align->get_results();

                    QUERY_FLAG_CHANGE(QUERY_ONT_INTERPRO_GO, !results->parsed_go.empty());
                    QUERY_FLAG_CHANGE(QUERY_ONT_INTERPRO_PATHWAY, !results->pathways.empty());

                    if (!results->parsed_go.empty())
                        QUERY_FLAG_SET(QUERY_ONE_GO);
                    if (!results->pathways.empty())
                        QUERY_FLAG_SET(QUERY_ONE_KEGG);

                    break;
                }
                default:
                    return;
            }
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

QueryAlignment* QuerySequence::AlignmentData::get_best_align_ptr(ExecuteStates state, uint16 software, std::string database) {
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

bool QuerySequence::AlignmentData::sort_descending_database::operator()(QueryAlignment *first,
                                                                        QueryAlignment *second) {
        first->set_compare_overall_alignment(false);
        second->set_compare_overall_alignment(false);
        return *first > *second;
}
