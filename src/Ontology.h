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

#ifndef ENTAP_ONTOLOGY_H
#define ENTAP_ONTOLOGY_H

#include "common.h"
#include <boost/program_options/variables_map.hpp>
#include "EntapConfig.h"
#include "EntapGlobals.h"
#include "QuerySequence.h"
#include "database/SQLDatabaseHelper.h"
#include "GraphingManager.h"
#include "ontology/AbstractOntology.h"
#include "QueryData.h"
#include "EntapModule.h"

class Ontology {

public:

    void execute();
    Ontology(std::string, EntapDataPtrs &);

private:

    const std::string ONTOLOGY_OUT_PATH     = "ontology/";
    const std::string FINAL_ANNOT_FILE      = "final_annotations_lvl";
    const std::string FINAL_ANNOT_FILE_CONTAM = "_contam";
    const std::string FINAL_ANNOT_FILE_NO_CONTAM = "_no_contam";
    const std::string ANNOT_FILE_EXT        = ".tsv";
    const uint16      FINAL_ALL_IND         = 0;
    const uint16      FINAL_CONTAM_IND      = 1;
    const uint16      FINAL_NO_CONTAM_IND   = 2;

    std::vector<std::string>        _interpro_databases;
    std::vector<uint16>             _go_levels;
    int                              _threads;
    std::vector<uint16>             _software_flags;
    bool                            _is_overwrite;
    bool                            _blastp;
    std::string                     _outpath;
    std::string                     _new_input;
    std::string                     _ontology_dir;
    std::string                     _eggnog_db_path;
    std::string                     _final_outpath_dir;
    std::vector<const std::string*> _HEADERS;
    GraphingManager                 *_pGraphingManager;
    QueryData                       *_QUERY_DATA;
    FileSystem                      *_pFileSystem;
    UserInput                       *_pUserInput;
    EntapDatabase                   *_pEntapDatabase;
    EntapDataPtrs                   _entap_data_ptrs;

    void print_eggnog(QUERY_MAP_T&);
    void init_headers();
    std::unique_ptr<EntapModule> spawn_object(uint16&);
};


#endif //ENTAP_ONTOLOGY_H
