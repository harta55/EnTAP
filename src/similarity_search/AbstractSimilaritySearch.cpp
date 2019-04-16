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

#include <regex>
#include "AbstractSimilaritySearch.h"

AbstractSimilaritySearch::AbstractSimilaritySearch(std::string &execution_stage_path, std::string &in_hits,
                                                   EntapDataPtrs &entap_data, std::string mod_name,
                                                   std::string &exe, vect_str_t &databases)
    :EntapModule(execution_stage_path, in_hits, entap_data, mod_name, exe) {

    _execution_state = SIMILARITY_SEARCH;

    // Get relevant user info for similarity searching
    _input_species    = _pUserInput->get_target_species_str();
    _qcoverage        = _pUserInput->get_user_input<fp32>(_pUserInput->INPUT_FLAG_QCOVERAGE);
    _tcoverage        = _pUserInput->get_user_input<fp32>(_pUserInput->INPUT_FLAG_TCOVERAGE);
    _e_val            = _pUserInput->get_user_input<fp64>(_pUserInput->INPUT_FLAG_E_VAL);
    _contaminants     = _pUserInput->get_contaminants();
    _uninformative_vect= _pUserInput->get_uninformative_vect();

    _database_paths = databases;

    // Get input species lineage information
    TaxEntry taxEntry = _pEntapDatabase->get_tax_entry(_input_species);
    _input_lineage    = taxEntry.lineage;

    // set blast string to use for file naming
    _blastp ? _blast_type = BLASTP_STR : _blast_type = BLASTX_STR;

    // create overall results dir
    _pFileSystem->delete_dir(_overall_results_dir);
    _pFileSystem->create_dir(_overall_results_dir);

}

// Returns database "shortname" from database full path
// This is just the filename, without any extension or path
std::string AbstractSimilaritySearch::get_database_shortname(std::string &full_path) {
    return _pFileSystem->get_filename(full_path, false);
}

std::string AbstractSimilaritySearch::get_database_output_path(std::string &database_name) {

    return PATHS(_blast_type + "_" + _transcript_shortname + "_" + get_database_shortname(database_name) + FileSystem::EXT_OUT,
            _mod_out_dir);
}

std::pair<bool,std::string> AbstractSimilaritySearch::is_contaminant(std::string lineage, vect_str_t &contams) {
    // species and tax database both lowercase
    if (contams.empty()) return std::pair<bool,std::string>(false,"");
    LOWERCASE(lineage);
    for (auto const &contaminant:contams) {
        if (lineage.find(contaminant) != std::string::npos){
            return std::pair<bool,std::string>(true,contaminant);
        }
    }
    return std::pair<bool,std::string>(false,"");
}

bool AbstractSimilaritySearch::is_informative(std::string title, vect_str_t &uninformative_vect) {
    LOWERCASE(title);
    for (std::string &item : uninformative_vect) { // Already lowercase
        if (title.find(item) != std::string::npos) return false;
    }
    return true;
}

std::string AbstractSimilaritySearch::get_species(std::string &title) {
    // TODO fix issue

    std::string species;


#ifdef USE_BOOST
    boost::smatch match;

    boost::regex ncbi_exp(_NCBI_REGEX);
    boost::regex uniprot_exp(_UNIPROT_REGEX);

    if (boost::regex_search(title,match,uniprot_exp)) {
        species = std::string(match[1].first, match[1].second);
    } else {
        if (boost::regex_search(title, match, ncbi_exp))
            species = std::string(match[1].first, match[1].second);
    }
#else   // Use std c++ libs
    std::smatch match;

    // Check if UniProt match
    if (std::regex_search(title, match, std::regex(UNIPROT_REGEX)) && match.size() > 1) {
        species = std::string(match[1].first, match[1].second);
    } else {
        // Not a UniProt match, check NCBI standard
        if (std::regex_search(title, match, std::regex(NCBI_REGEX)) && match.size() > 1) {
            species = std::string(match[1].first, match[1].second);
        }
    }
#endif

    // Double bracket fix
    if (species[0] == '[') species = species.substr(1);
    if (species[species.length()-1] == ']') species = species.substr(0,species.length()-1);

    return species;
}

bool AbstractSimilaritySearch::is_uniprot_entry(std::string &sseqid, UniprotEntry &entry) {
    std::string accession;

    // sseqid - sp|Q9FJZ9|PER72_ARATH
    if (_pEntapDatabase == nullptr || sseqid.empty()) return false;

    accession = sseqid.substr(sseqid.rfind('|',sseqid.length())+1);     // Q9FJZ9
    entry = _pEntapDatabase->get_uniprot_entry(accession);
    return !entry.is_empty();
}