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

//*********************** Includes *****************************
#include "SimilaritySearch.h"
#include "similarity_search/ModDiamond.h"
//**************************************************************

/**
 * ======================================================================
 * Function SimilaritySearch(std::vector<std::string> &databases, std::string input,
                           int threads, std::string out, boost::program_options::variables_map &user_flags,
                           GraphingManager *graphingManager)
 *
 * Description          - SimilaritySearch object constructor
 *                      - Responsible for initiating SimSearch member variables,
 *                        parsing user input for relevant information used within
 *                        this module
 *
 * Notes                - None
 *
 * @param databases     - List of user selected databases (parsed in  execute namespace)
 * @param input         - Path to user transcriptome (final version post filtering)
 * @param threads       - Thread count
 * @param out           - EnTAP out directory
 * @param user_flags    - Boost parsed user input
 * @param graphingManager- Pointer to graphing manager
 *
 * @return              - SimilaritySearch instance
 * ======================================================================
 */
SimilaritySearch::SimilaritySearch(vect_str_t &database_paths, std::string &input, EntapDataPtrs &entap_data) {
    FS_dprint("Spawn Object - SimilaritySearch");
    std::string uninform_path;

    mpUserInput     = entap_data.mpUserInput;
    mpFileSystem    = entap_data.mpFileSystem;
    mpEntapData    = &entap_data;
    mInputFastaPath     = input;
    mDatabasePaths = database_paths;

    // Set sim search paths/directories
    mSimSearchDir  = PATHS(mpFileSystem->get_root_path(), SIM_SEARCH_DIR);
    mOverwrite = mpUserInput->has_input(INPUT_FLAG_OVERWRITE);

    if (mOverwrite) {
        mpFileSystem->delete_dir(mSimSearchDir);
    }
    mpFileSystem->create_dir(mSimSearchDir);
    mSoftwareFlag = SIM_DIAMOND;
}

void SimilaritySearch::execute() {
    EntapModule::ModVerifyData verifyData;
    std::unique_ptr<AbstractSimilaritySearch> ptr;

    try {
        ptr = spawn_object();
        verifyData = ptr->verify_files();
        if (!verifyData.files_exist) {
            ptr->execute();
        }
        ptr->parse();
        ptr->set_success_flags();

        ptr.reset();
    } catch (const ExceptionHandler &e) {
        ptr.reset();
        throw e;
    }
}

std::unique_ptr<AbstractSimilaritySearch> SimilaritySearch::spawn_object() {

    switch (mSoftwareFlag) {
        // Fall through to default case (only 1 supported Sim Search module)
        case SIM_DIAMOND:
        default:
            return std::unique_ptr<AbstractSimilaritySearch>(new ModDiamond(
                    mSimSearchDir, mInputFastaPath, *mpEntapData, mDatabasePaths
                    ));
    }
}

