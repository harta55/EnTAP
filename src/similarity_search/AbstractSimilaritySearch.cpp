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

#include <regex>
#include "AbstractSimilaritySearch.h"
#include "../QueryData.h"

/**
 * ======================================================================
 * Function AbstractSimilaritySearch::AbstractSimilaritySearch(std::string &execution_stage_path, std::string &in_hits,
 *                                     EntapDataPtrs &entap_data, std::string mod_name,
 *                                     std::string &exe, vect_str_t &databases)
 *
 * Description           - Type of EntapModule that handles Similarity Searching
 *
 * Notes                 - Inherits from EntapModule
 *
 * @param execution_stage_path - Absolute directory path to execution stage (Similarity Search)
 * @param in_hits              - Absolute path to FASTA transcriptome
 * @param entap_data           - Structure of necessary pointers for execution
 * @param mod_name             - Name of module (will be used for directory naming)
 * @param exe                  - Module execution method (DIAMOND exe file)
 * @param databases            - Vector of databases to search against
 *
 * @return               - AbstractSimilaritySearch object
 *
 * =====================================================================
 */
AbstractSimilaritySearch::AbstractSimilaritySearch(std::string &execution_stage_path, std::string &in_hits,
                                                   EntapDataPtrs &entap_data, std::string mod_name,
                                                   std::vector<ENTAP_HEADERS> &module_headers,
                                                   vect_str_t &databases)
    :EntapModule(execution_stage_path, in_hits, entap_data, mod_name, module_headers) {

    mExecutionState = SIMILARITY_SEARCH;

    // Get relevant user info for similarity searching
    mInputSpecies    = mpUserInput->get_target_species_str();
    mQCoverage        = mpUserInput->get_user_input<ent_input_fp_t >(INPUT_FLAG_QCOVERAGE);
    mTCoverage        = mpUserInput->get_user_input<ent_input_fp_t >(INPUT_FLAG_TCOVERAGE);
    mEVal            = mpUserInput->get_user_input<ent_input_fp_t>(INPUT_FLAG_E_VALUE);
    mContaminateTaxons     = mpUserInput->get_contaminants();
    mUninformativeTags= mpUserInput->get_uninformative_vect();

    mDatabasePaths = databases;

    // Get input species lineage information
    TaxEntry taxEntry = mpEntapDatabase->get_tax_entry(mInputSpecies);
    mInputLineage    = taxEntry.lineage;

    // set blast string to use for file naming
    mBlastp ? mBlastType = BLASTP_STR : mBlastType = BLASTX_STR;

    // create overall results dir
    mpFileSystem->delete_dir(mOverallResultsDir);
    mpFileSystem->create_dir(mOverallResultsDir);

}

/**
 * ======================================================================
 * Function std::string AbstractSimilaritySearch::get_database_shortname
 *                      (std::string &full_path)
 *
 * Description           - Returns database 'shortname' from full database path
 *                         (filename without extension)
 *
 * Notes                 - None
 *
 * @param full_path      - Absolute path to database
 *
 * @return               - Database shortname (filename without extension)
 *
 * =====================================================================
 */
std::string AbstractSimilaritySearch::get_database_shortname(std::string &full_path) {
    return mpFileSystem->get_filename(full_path, false);
}

/**
 * ======================================================================
 * Function std::string AbstractSimilaritySearch::get_database_output_path
 *                                              (std::string &database_name)
 *
 * Description           - Gets the desired outpath from DIAMOND/blast execution
 *
 * Notes                 - None
 *
 * @param database_name  - Name of database to execute DIAMOND/BLAST against
 *
 * @return               - Absolute path to desired output path post-sim search
 *
 * =====================================================================
 */
std::string AbstractSimilaritySearch::get_database_output_path(std::string &database_name) {

    return PATHS(mModOutDir,mBlastType + "_" + mTranscriptomeShortname + "_" + get_database_shortname(database_name) + FileSystem::EXT_OUT);
}

/**
 * ======================================================================
 * Function std::pair<bool,std::string> AbstractSimilaritySearch::is_contaminant
 *                                      (std::string lineage, vect_str_t &contams)
 *
 * Description           - Determine if specified lineage is a contaminant based
 *                         on input contaminant taxons
 *
 * Notes                 - None
 *
 * @param lineage        - Target lineage to determine contaminant status
 * @param contams        - vector of contaminant taxons
 *
 * @return               - BOOL (contaminate TRUE/FALSE) and STRING (matching taxon)
 *
 * =====================================================================
 */
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

/**
 * ======================================================================
 * Function bool AbstractSimilaritySearch::is_informative(std::string title,
 *                                          vect_str_t &uninformative_vect)
 *
 * Description           - Determine if specified subject title is 'informative'
 *                         based upon informative tags
 *
 * Notes                 - None
 *
 * @param title              - Target title to determine informative status
 * @param uninformative_vect - vector of 'informative' tags
 *
 * @return               - TRUE/FALSE informative
 *
 * =====================================================================
 */
bool AbstractSimilaritySearch::is_informative(std::string title, vect_str_t &uninformative_vect) {
    LOWERCASE(title);
    for (std::string &item : uninformative_vect) { // Already lowercase
        if (title.find(item) != std::string::npos) return false;
    }
    return true;
}

/**
 * ======================================================================
 * Function std::string AbstractSimilaritySearch::get_species(std::string &title)
 *
 * Description           - Pull species from NCBI/UniProt or other title format
 *                       - BOOST/regex support
 *
 * Notes                 - None
 *
 * @param title          - Target title to parse for species
 *
 * @return               - Species parsed from title input
 *
 * =====================================================================
 */
std::string AbstractSimilaritySearch::get_species(std::string &title) {
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
        // GCC 4.8.x does NOT fully support std regex and have bugs, full...bugless.. starts at 5.x
#if 0
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
    uint64 ind1, ind2;
    ind1 = title.find("OS=");

    // if UniProt type format for species EX: OS=Homo sapiens
    if (ind1 != std::string::npos) {
        // Yes,
        ind2 = title.find('=', ind1 + 3);
        if (ind2 != std::string::npos && (ind2 - ind1) > 6) {
            species = title.substr(ind1 + 3, ind2 - ind1 - 6);
        }

    } else {
        // No, possibly NCBI format. EX: [homo sapiens]
        ind1 = title.find_last_of('[');
        ind2 = title.find_last_of(']');
        if (ind1 != std::string::npos && ind2 != std::string::npos) {
            species = title.substr(ind1 + 1, (ind2 - ind1 - 1));
        }
    }
#endif

    // Double bracket fix
    if (species[0] == '[') species = species.substr(1);
    if (species[species.length()-1] == ']') species = species.substr(0,species.length()-1);

    return species;
}

void AbstractSimilaritySearch::set_success_flags() {
    mpQueryData->set_is_success_sim_search(true);
}
