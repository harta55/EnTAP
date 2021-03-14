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

#ifndef ENTAP_SIMILARITYSEARCH_H
#define ENTAP_SIMILARITYSEARCH_H

//*********************** Includes *****************************
#include "common.h"
#include "QueryData.h"
#include "database/EntapDatabase.h"
#include "FileSystem.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
#include "UserInput.h"
#include "database/EntapDatabase.h"
#include "similarity_search/AbstractSimilaritySearch.h"

//**************************************************************


class SimilaritySearch {

public:

    //******************** Public Prototype Functions *********************
    SimilaritySearch(vect_str_t &database_paths, std::string &input, EntapDataPtrs &entapDataPtrs);
    void execute();


    //*********************************************************************

private:

    vect_str_t                      mDatabasePaths;     // Absolute paths to databases input from user
    std::string                     mInputFastaPath;    // Absolute path to input FASTA file from user
    std::string                     mSimSearchDir;      // Absolute path to root Similarity Search directory
    SIMILARITY_SOFTWARE             mSoftwareFlag;      // Selected sim search module
    bool                            mOverwrite;         // TRUE if user would like to remove previous data
    FileSystem                      *mpFileSystem;      // Pointer to EnTAP File System
    UserInput                       *mpUserInput;       // Pointer to User Input
    EntapDataPtrs                   *mpEntapData;       // Pointer to all data required for execution

    std::string SIM_SEARCH_DIR = "similarity_search";

    std::unique_ptr<AbstractSimilaritySearch> spawn_object();
};


#endif //ENTAP_SIMILARITYSEARCH_H
