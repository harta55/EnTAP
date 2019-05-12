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

    vect_str_t                      _database_paths;
    std::string                     _diamond_exe;
    std::string                     _outpath;
    std::string                     _input_path;
    std::string                     _sim_search_dir;
    SIMILARITY_SOFTWARE             _software_flag;
    bool                            _overwrite;
    FileSystem                      *_pFileSystem;
    UserInput                       *_pUserInput;
    EntapDataPtrs                   *_pEntap_data;

    std::string SIM_SEARCH_DIR = "similarity_search";

    std::unique_ptr<AbstractSimilaritySearch> spawn_object();
};


#endif //ENTAP_SIMILARITYSEARCH_H
