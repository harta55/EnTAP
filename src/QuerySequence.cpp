/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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
#include "common.h"
#include "QueryAlignment.h"

/**
 * ======================================================================
 * Function uint64 QuerySequence::get_sequence_length() const
 *
 * Description           - Returns nucleotide sequence length
 *
 * Notes                 - None
 *
 * @return               - Sequence length (nucleotide base pairs)
 *
 * =====================================================================
 */
uint64 QuerySequence::get_sequence_length() const {
    return mSequenceLength;
}

/**
 * ======================================================================
 * Function QuerySequence::QuerySequence()
 *
 * Description           - Initializes QuerySequence object
 *
 * Notes                 - Constructor
 *
 * @return               - Query Sequence object
 *
 * =====================================================================
 */
QuerySequence::QuerySequence() {
    init_sequence();
}

/**
 * ======================================================================
 * Function QuerySequence::QuerySequence(bool is_protein, std::string seq,
 *                                       std::string seqid)
 *
 * Description          - Initializes QuerySequence object with sequence, ID, and length
 *                      - Calculates sequence length and may set QUERY_IS_PROTEIN
 *                        flag
 *
 * Notes                - Constructor
 *
 * @param is_protein    - TRUE/FALSE if protein/nucleotide sequence
 * @param seq           - Sequence to set for QuerySequence
 * @param seqid         - Sequence ID
 *
 *
 * @return              - QuerySequence object
 *
 * =====================================================================
 */
QuerySequence::QuerySequence(bool is_protein, std::string seq, std::string seqid) {
    init_sequence();
    trim_sequence(seq);

    mSequenceID = seqid;

    if (is_protein) {
        QUERY_FLAG_SET(QUERY_IS_PROTEIN);
        mSequenceProtein = seq;
    } else {
        QUERY_FLAG_CLEAR(QUERY_IS_PROTEIN);
        QUERY_FLAG_SET(QUERY_IS_NUCLEOTIDE);
        mSequenceNucleo = seq;
    }

    mSequenceLength = calc_seq_length(seq,is_protein);
    set_header_data();
}

/**
 * ======================================================================
 * Function QuerySequence::~QuerySequence()
 *
 * Description           - Cleans up allocated memory (all alignment data)
 *
 * Notes                 - Destructor
 *
 * @return               - None
 *
 * =====================================================================
 */
QuerySequence::~QuerySequence() {
    // Clear alignment data
    delete mAlignmentData;
}

/**
 * ======================================================================
 * Function const std::string &QuerySequence::get_sequence_p() const
 *
 * Description          - Returns protein sequence (amino acid)
 *
 * Notes                - None
 *
 * @return              - Protein sequence
 *
 * =====================================================================
 */
const std::string &QuerySequence::get_sequence_p() const {
    return mSequenceProtein;
}

/**
 * ======================================================================
 * Function void QuerySequence::set_sequence_p(std::string &seq)
 *
 * Description          - Sets protein (amino acid) sequence and re-calculates
 *                        sequence length (bp)
 *                      - Sets QUERY_IS_PROTEIN flag
 *
 * Notes                - None
 *
 * @param seq           - Protein (amino acid) sequence
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::set_sequence_p(std::string &seq) {
    QUERY_FLAG_SET(QUERY_IS_PROTEIN);
    mSequenceLength = calc_seq_length(seq,true);
    if (!seq.empty() && seq[seq.length()-1] == '\n') {
        seq.pop_back();
    }
    mSequenceProtein = seq;
}

/**
 * ======================================================================
 * Function const std::string &QuerySequence::get_sequence_n() const
 *
 * Description          - Returns nucleotide sequence
 *
 * Notes                - None
 *
 * @return              - Nucleotide sequence
 *
 * =====================================================================
 */
const std::string &QuerySequence::get_sequence_n() const {
    return mSequenceNucleo;
}

/**
 * ======================================================================
 * Function void QuerySequence::set_sequence_n(const std::string &_sequence_n)
 *
 * Description          - Sets nucleotide sequence
 *
 * Notes                - None
 *
 * @param sequence_n    - Nucleotide sequence
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::set_sequence_n(std::string &sequence_n) {
    QUERY_FLAG_SET(QUERY_IS_NUCLEOTIDE);
    if (!sequence_n.empty() && sequence_n[sequence_n.length()-1] == '\n') {
        sequence_n.pop_back();
    }
    QuerySequence::mSequenceNucleo = sequence_n;
}

/**
 * ======================================================================
 * Function void QuerySequence::trim_sequence(std::string &sequence)
 *
 * Description          - Removes trailing newline character
 *
 * Notes                - None
 *
 * @param sequence      - Sequence to be trimmed
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::trim_sequence(std::string &sequence) {
    if (!sequence.empty() && sequence[sequence.length()-1] == '\n') {
        sequence.pop_back();
    }
}

/**
 * ======================================================================
 * Function unsigned long QuerySequence::calc_seq_length(std::string &seq,
 *                                                      bool protein)
 *
 * Description          - Calculates sequence length in bp from nucleotide
 *                        or protein input sequence
 *
 * Notes                - None
 *
 * @param protein       - TRUE/FALSE if protein/nucleotide sequence
 * @param seq           - Sequence to calculate bp length for
 *
 *
 * @return              - Sequence bp length
 *
 * =====================================================================
 */
uint64 QuerySequence::calc_seq_length(std::string &seq,bool protein) {
    std::string sub = seq.substr(seq.find('n')+1);
    long line_chars = std::count(sub.begin(),sub.end(),'\n');
    uint64 seq_len = sub.length() - line_chars;
    if (protein) seq_len *= NUCLEO_PER_AMINO;
    return seq_len;
}

/**
 * ======================================================================
 * Function const std::string &QuerySequence::getFrame() const
 *
 * Description          - Returns frame type of sequence (internal,
 *                        complete, partial...etc)
 *
 * Notes                - None
 *
 *
 * @return              - Type of frame for sequence (partial, internal,etc)
 *
 * =====================================================================
 */
const std::string &QuerySequence::getFrame() const {
    return mFrameType;
}

/**
 * ======================================================================
 * Function void QuerySequence::setFrame(const std::string &frame)
 *
 * Description          - Sets frame type (Defined in ModGeneMarkST.h)
 *
 * Notes                - None
 *
 * @param frame         - Frame type (defined in ModGeneMarkST.h)
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::setFrame(const std::string &frame) {
    QuerySequence::mFrameType = frame;
    set_header_data();
}

#ifdef EGGNOG_MAPPER
void QuerySequence::set_eggnog_results(const EggnogResults &eggnogResults) {
    memcpy(&this->mEggnogResults, &eggnogResults, sizeof(eggnogResults));
    QUERY_FLAG_SET(QUERY_EGGNOG_HIT);
    QUERY_FLAG_SET(QUERY_FAMILY_ASSIGNED);
}
#endif

/**
 * ======================================================================
 * Function void QuerySequence::init_sequence()
 *
 * Description          - INITs sequence data
 *                      - Sets QUERY_FRAME/EXPRESSION_KEPT flags
 *                      - Initializes header mappings
 *
 * Notes                - None
 *
 * @param frame         - Frame type (defined in ModGeneMarkST.h)
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::init_sequence() {
    mSequenceLength = 0;
    mFPKM = 0.0;
    mTPM = 0.0;
    mEffectiveLength = 0.0;
    mFrameScore = 0.0;

    mAlignmentData = new AlignmentData(this);
#ifdef EGGNOG_MAPPER
    mEggnogResults = EggnogResults();
#endif

    mFrameType = "";
    mSequenceProtein = "";
    mSequenceNucleo = "";

    mQueryFlags = 0;
    QUERY_FLAG_SET(QUERY_FRAME_KEPT);
    QUERY_FLAG_SET(QUERY_EXPRESSION_KEPT);
    set_header_data();
}

/**
 * ======================================================================
 * Function const std::string &QuerySequence::get_sequence() const
 *
 * Description          - Returns nucleotide sequence (if exists) or protein
 *
 * Notes                - Used during Expression Analysis only
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
const std::string &QuerySequence::get_sequence() const {
    if (mSequenceNucleo.empty()) {
        return mSequenceProtein;
    }
    else {
        return mSequenceNucleo;
    }
}

/**
 * ======================================================================
 * Function void QuerySequence::get_header_data(std::string &data, ENTAP_HEADERS header,
 *                                              uint8 lvl)
 *
 * Description          - Returns data corresponding to the requested header from
 *                        mHeaderInfo
 *                      - May access alignment data depending on header
 *                      - Normalizes GO terms to input level
 *
 * Notes                - None
 *
 * @param data          - Pointer to requested header data, will be modified
 * @param header        - Requested header
 * @param lvl           - GO term level to normalize to
 *
 * @return              - None (modifies data ptr)
 *
 * =====================================================================
 */
void QuerySequence::get_header_data(std::string &data, ENTAP_HEADERS header, uint8 lvl) {

    QueryAlignment *align_ptr;      // Pointer to alignment data

    align_ptr = nullptr;
    data = "";

    // Need to re-access alignment data to normalize to the GO level
    // Otherwise, access member header data
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
            data = mHeaderInfo[header];
            break;
    }

    if (align_ptr != nullptr) {
        align_ptr->get_header_data(header, data, lvl);
    }
}

/**
 * ======================================================================
 * Function void QuerySequence::set_header_data()
 *
 * Description          - Updates the header data in member mHeaderInfo variable
 *                      - Accesses "best hit" for alignment data and updates header
 *                        info
 *
 * Notes                - None
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::set_header_data() {
    QueryAlignment *align_ptr = nullptr;

    // General data
    mHeaderInfo[ENTAP_HEADER_QUERY] = this->mSequenceID;

    // Frame Selection data
    mHeaderInfo[ENTAP_HEADER_FRAME] = this->mFrameType;

    // Expression Filtering data
    mHeaderInfo[ENTAP_HEADER_EXP_FPKM] = float_to_string(this->mFPKM);
    mHeaderInfo[ENTAP_HEADER_EXP_TPM] = float_to_string(this->mTPM);
    mHeaderInfo[ENTAP_HEADER_EXP_E_LENGTH] = float_to_string(this->mEffectiveLength);

    // Similarity Search data
    align_ptr = this->mAlignmentData->get_best_align_ptr(SIMILARITY_SEARCH, SIM_DIAMOND, "");
    if (align_ptr != nullptr) {
        align_ptr->get_all_header_data(mHeaderInfo);
    }

    // Ontology EggNOG data
    align_ptr = this->mAlignmentData->get_best_align_ptr(GENE_ONTOLOGY, ONT_EGGNOG_DMND, "");
    if (align_ptr != nullptr) {
        align_ptr->get_all_header_data(mHeaderInfo);
    }

    // Ontology InterProScan data
    align_ptr = this->mAlignmentData->get_best_align_ptr(GENE_ONTOLOGY, ONT_INTERPRO_SCAN, "");
    if (align_ptr != nullptr) {
        align_ptr->get_all_header_data(mHeaderInfo);
    }

    // Ontology BUSCO data
    align_ptr = this->mAlignmentData->get_best_align_ptr(GENE_ONTOLOGY, ONT_BUSCO, "");
    if (align_ptr != nullptr) {
        align_ptr->get_all_header_data(mHeaderInfo);
    }
}

/**
 * ======================================================================
 * Function void QuerySequence::set_fpkm(float fpkm)
 *
 * Description          - Sets the FPKM value from Expression Filtering
 *
 * Notes                - None
 *
 * @param fpkm          - FPKM value to set
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::set_fpkm(fp32 fpkm) {
    mFPKM = fpkm;
    set_header_data();  // Reset header data
}

/**
 * ======================================================================
 * Function bool QuerySequence::isContaminant()
 *
 * Description          - Returns contaminant status (from Similarity Search)
 *
 * Notes                - None
 *
 *
 * @return              - TRUE/FALSE whether sequence was flagged as a
 *                        contaminant
 *
 * =====================================================================
 */
bool QuerySequence::is_contaminant() {
    return this->QUERY_FLAG_GET(QUERY_CONTAMINANT);
}

/**
 * ======================================================================
 * Function bool QuerySequence::is_kept()
 *
 * Description          - Returns status if the sequence has not been removed
 *                        through Expression Filtering or Frame Selection
 *
 * Notes                - None
 *
 *
 * @return              - TRUE/FALSE whether sequence is kept or removed
 *                        from analysis
 *
 * =====================================================================
 */
bool QuerySequence::is_kept() {
    return QUERY_FLAG_GET(QUERY_EXPRESSION_KEPT) && QUERY_FLAG_GET(QUERY_FRAME_KEPT);
}


bool QuerySequence::QUERY_FLAG_GET(QUERY_FLAGS flag) {
    return (mQueryFlags & flag) != 0;
}

void QuerySequence::QUERY_FLAG_SET(QUERY_FLAGS flag) {
    mQueryFlags |= flag;
}

void QuerySequence::QUERY_FLAG_CLEAR(QUERY_FLAGS flag) {
    mQueryFlags &= ~flag;
}

void QuerySequence::QUERY_FLAG_CHANGE(QuerySequence::QUERY_FLAGS flag, bool val) {
    if (val) {
        QUERY_FLAG_SET(flag);
    } else {
        QUERY_FLAG_CLEAR(flag);
    }
}

/**
 * ======================================================================
 * Function void QuerySequence::add_alignment(ExecuteStates state, uint16 software,
 *                                      EggnogResults &results, std::string& database)
 *
 * Description          - Adds EggNOG alignment to AlignmentData and updates
 *                        pertinent data/best hits (flags, EggnogResults...)
 *
 * Notes                - None
 *
 * @param state         - State that alignment was created in
 * @param software      - Software that alignment was created using (DIAMOND, EggNOG...)
 * @param results       - EggNOG data
 * @param database      - Absolute path to database associated with alignment
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::add_alignment(ExecuteStates state, uint16 software, EggnogResults &results, std::string& database) {
    QUERY_FLAG_SET(QUERY_EGGNOG_HIT);
    QUERY_FLAG_SET(QUERY_FAMILY_ASSIGNED);
    QueryAlignment *new_alignment = new EggnogDmndAlignment(state, software, database, this, results);
    mAlignmentData->update_best_hit(new_alignment);
}

/**
 * ======================================================================
 * Function void QuerySequence::add_alignment(ExecuteStates state, uint16 software,
 *                                      SimSearchResults &results, std::string& database)
 *
 * Description          - Adds Similarity Search alignment to AlignmentData and updates
 *                        pertinent data/best hits (flags, EggnogResults...)
 *
 * Notes                - None
 *
 * @param state         - State that alignment was created in
 * @param software      - Software that alignment was created using (DIAMOND, EggNOG...)
 * @param results       - Similarity Search data
 * @param database      - Absolute path to database associated with alignment
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::add_alignment(ExecuteStates state, uint16 software, SimSearchResults &results, std::string& database,std::string lineage) {
    QUERY_FLAG_SET(QUERY_BLAST_HIT);
    QueryAlignment *new_alignment = new SimSearchAlignment(state, software, database, this, results, lineage);
    mAlignmentData->update_best_hit(new_alignment);
}

/**
 * ======================================================================
 * Function void QuerySequence::add_alignment(ExecuteStates state, uint16 software,
 *                                      InterProResults &results, std::string& database)
 *
 * Description          - Adds InterPro alignment to AlignmentData and updates
 *                        pertinent data/best hits (flags, EggnogResults...)
 *
 * Notes                - None
 *
 * @param state         - State that alignment was created in
 * @param software      - Software that alignment was created using (DIAMOND, EggNOG...)
 * @param results       - InterPro data
 * @param database      - Absolute path to database associated with alignment
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::add_alignment(ExecuteStates state, uint16 software, InterProResults &results,
                                  std::string &database) {
    QUERY_FLAG_SET(QUERY_INTERPRO);
    QueryAlignment *new_alignmet = new InterproAlignment(state, software, database, this, results);
    mAlignmentData->update_best_hit(new_alignmet);
}


/**
 * ======================================================================
 * Function void QuerySequence::add_alignment(ExecuteStates state, uint16 software,
 *                                      BuscoResults &results, std::string& database)
 *
 * Description          - Adds BuscoAlignment alignment to AlignmentData and updates
 *                        pertinent data/best hits
 *
 * Notes                - None
 *
 * @param state         - State that alignment was created in
 * @param software      - Software that alignment was created using (DIAMOND, EggNOG...)
 * @param results       - BUSCO data
 * @param database      - Absolute path to database associated with alignment
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::add_alignment(ExecuteStates state, uint16 software, QuerySequence::BuscoResults &results,
                                  std::string &database) {

    QUERY_FLAG_SET(QUERY_ONT_BUSCO);
    QueryAlignment *new_alignmet = new BuscoAlignment(state, software, database, this, results);
    mAlignmentData->update_best_hit(new_alignmet);
}

//**********************************************************************
//**********************************************************************
//                 AlignmentData Nested Structure
//**********************************************************************
//**********************************************************************

/**
 * ======================================================================
 * Function QuerySequence::AlignmentData::AlignmentData(QuerySequence *sequence)
 *
 * Description          - Initializes alignment data class
 *
 * Notes                - Constructor
 *
 * @param sequence      - Pointer to parent sequence
 *
 * @return              - AlignmentData object
 *
 * =====================================================================
 */
QuerySequence::AlignmentData::AlignmentData(QuerySequence *sequence){
    querySequence = sequence;
}

/**
 * ======================================================================
 * Function QuerySequence::AlignmentData::~AlignmentData()
 *
 * Description          - Cleans up AlignmentData memory by accessing
 *                        all QueryAlignments
 *
 * Notes                - Destructor
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
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

/**
 * ======================================================================
 * Function void QuerySequence::AlignmentData::update_best_hit(ExecuteStates state,
 *          uint16 software, std::string &database, QueryAlignment* new_alignment)
 *
 * Description          - Updates all alignment data including best hit for
 *                        each database, overall software best hit across
 *                        every database, calls Query flag update
 *
 * Notes                - None
 *
 * @param state         - Execution state the algnment belongs to
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void QuerySequence::AlignmentData::update_best_hit(QueryAlignment* new_alignment) {

    ALIGNMENT_DATA_T* alignment_arr = get_software_ptr(new_alignment->getMExecutionState(), new_alignment->getMSoftwareModule());

    // Did we hit against this database yet
    if (!hit_database(new_alignment->getMExecutionState(), new_alignment->getMSoftwareModule(),
                      new_alignment->getMDatabasePath())) {
        // No, create new vector for that database and add as best hit for database
        align_database_hits_t vect = {new_alignment};
        alignment_arr->emplace(new_alignment->getMDatabasePath(), vect);
    } else {
        // Yes, add alignment to list then update
        alignment_arr->at(new_alignment->getMDatabasePath()).push_back(new_alignment);
    }

    // Always will have hit this database, get alignments
    align_database_hits_t *database_data = &alignment_arr->at(new_alignment->getMDatabasePath());

    // Sort alignments for that database (0 index is best hit)
    std::sort(database_data->begin(), database_data->end(), sort_descending_database());

    // See if this alignment is better than the overall alignment
    new_alignment->set_compare_overall_alignment(true);
    QueryAlignment* best_alignment = get_best_align_ptr(new_alignment->getMExecutionState(), new_alignment->getMSoftwareModule(), "");

    if (best_alignment == nullptr) {
        set_best_alignment(new_alignment);
    } else {
        best_alignment->set_compare_overall_alignment(true);
        if (*new_alignment > *best_alignment) {
            set_best_alignment(new_alignment);
        }
    }

    // Update any overall flags that may have changed with best hit changes
    querySequence->update_query_flags(new_alignment->getMExecutionState(), new_alignment->getMSoftwareModule());
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

                case ONT_BUSCO: {
                    // Nothing special to do here
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
    return this->mAlignmentData->get_database_ptr(state, software, database);
}

std::string QuerySequence::format_go_info(go_format_t &go_list, uint8 lvl) {
    std::stringstream out;
    std::string temp;

    for (GoEntry const &val : go_list)  {
        // Only return unknown levels when asking for '0' GO level
        // CLEANUP!
        if (lvl == 0 || (val.level_int >= lvl && val.level_int != GoEntry::UNKNOWN_LVL)) {
            temp = val.go_id + "(L=" + std::to_string(val.level_int) + ")";
            out << temp << ",";
        }
    }
    return out.str();
}

bool QuerySequence::hit_database(ExecuteStates state, uint16 software, std::string database) {
    return mAlignmentData->hit_database(state, software, database);
}

void QuerySequence::set_blasted() {
    QUERY_FLAG_SET(QUERY_BLASTED);
}

const std::string &QuerySequence::getMSequenceID() const {
    return mSequenceID;
}

bool QuerySequence::is_protein() {
    return QUERY_FLAG_GET(QUERY_IS_PROTEIN);
}

bool QuerySequence::is_kept_expression() {
    return QUERY_FLAG_GET(QUERY_EXPRESSION_KEPT);
}

void QuerySequence::setMTPM(fp64 mTPM) {
    QuerySequence::mTPM = mTPM;
}

bool QuerySequence::QUERY_FLAG_CONTAINS(uint32 flags) {
    return ((flags & mQueryFlags) != 0);
}

uint32 QuerySequence::getMQueryFlags() const {
    return mQueryFlags;
}

bool QuerySequence::is_nucleotide() {
    return QUERY_FLAG_GET(QUERY_IS_NUCLEOTIDE);
}

void QuerySequence::setMEffectiveLength(fp32 mEffectiveLength) {
    QuerySequence::mEffectiveLength = mEffectiveLength;
}

go_format_t QuerySequence::get_go_terms() {
    go_format_t  ret = go_format_t();
    go_format_t  align_data;

    // Pull EggNOG GO Terms
    auto *egg_alignment = get_best_hit_alignment<EggnogDmndAlignment>(GENE_ONTOLOGY, ONT_EGGNOG_DMND, "");
    if (egg_alignment != nullptr) {
        if (!egg_alignment->get_go_data().empty()) {
            ret = egg_alignment->get_go_data();
        }
    }

    // Pull UniProt GO Terms from Similarity Searching (overall best hit)
    auto *sim_alignment = get_best_hit_alignment<SimSearchAlignment>(SIMILARITY_SEARCH, SIM_DIAMOND, "");
    if (sim_alignment != nullptr) {
        align_data = sim_alignment->get_go_data();
        if (!align_data.empty()) {
            std::merge(ret.begin(), ret.end(), align_data.begin(), align_data.end(),
                       std::inserter(ret, ret.begin()));
        }
    }

    // Pull InterPro GO Terms
    auto *inter_alignment = get_best_hit_alignment<InterproAlignment>(GENE_ONTOLOGY, ONT_INTERPRO_SCAN, "");
    if (inter_alignment != nullptr) {
        align_data = inter_alignment->get_go_data();
        if (!align_data.empty()) {
            std::merge(ret.begin(), ret.end(), align_data.begin(), align_data.end(),
                       std::inserter(ret, ret.begin()));
        }
    }

    return ret;
}

fp32 QuerySequence::getMEffectiveLength() const {
    return mEffectiveLength;
}

bool QuerySequence::contains_go_level(int16 level) {
    go_format_t go_terms = get_go_terms();
    if (!go_terms.empty()) {
        if (level == 0) return true;
        for (GoEntry const &entry : go_terms) {
            if (entry.level_int == level) return true;
        }
    } else {
        return false;
    }

    return false;
}

fp32 QuerySequence::getMFrameScore() const {
    return mFrameScore;
}

void QuerySequence::setMFrameScore(fp32 mFrameScore) {
    QuerySequence::mFrameScore = mFrameScore;
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

QuerySequence::AlignmentData::ALIGNMENT_DATA_T* QuerySequence::AlignmentData::get_software_ptr(ExecuteStates state, uint16 software) {
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
QuerySequence::AlignmentData::set_best_alignment(QueryAlignment *new_alignment) {
    overall_alignment[new_alignment->getMExecutionState()][new_alignment->getMSoftwareModule()] = new_alignment;
}

bool QuerySequence::AlignmentData::sort_descending_database::operator()(QueryAlignment *first,
                                                                        QueryAlignment *second) {
        first->set_compare_overall_alignment(false);
        second->set_compare_overall_alignment(false);
        return *first > *second;
}