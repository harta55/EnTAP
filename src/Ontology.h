/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
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

#include <iostream>
#include <map>
#include <boost/program_options/variables_map.hpp>
#include "EntapConfig.h"
#include "EntapGlobals.h"
#include "QuerySequence.h"
#include "DatabaseHelper.h"
#include "GraphingManager.h"
#include "ontology/AbstractOntology.h"

class QuerySequence;
class AbstractOntology;
class QueryData;

class Ontology {
    typedef std::map<std::string, QuerySequence> query_map_struct;


public:

    void execute(std::string,std::string);
    Ontology(int,std::string,std::string, boost::program_options::variables_map &,GraphingManager*,
             QueryData*, bool);

private:

    enum OntologySoftware {
        EGGNOG,
        INTERPRO
    };

    OntologySoftware SOFTWARE;
    const std::string ONTOLOGY_OUT_PATH     = "ontology/";
    const std::string PROCESSED_OUT_DIR     = "processed/";
    const std::string FIGURE_DIR            = "figures/";

    std::vector<short>              _go_levels;
    int                             _threads;
    short                           _software_flag;
    bool                            _is_overwrite;
    bool                            _blastp;
    std::string                     _ontology_exe;
    std::string                     _outpath;
    std::string                     _new_input;
    std::string                     _input_no_hits;
    std::string                     _ontology_dir;
    std::string                     _processed_dir;
    std::string                     _figure_dir;
    std::string                     _eggnog_db_path;
    std::vector<const std::string*> _HEADERS;
    std::vector<std::string>        _interpro_databases;
    GraphingManager                 *_graphingManager;
    QueryData                       *_QUERY_DATA;

    void print_eggnog(query_map_struct&);
    void print_interpro(query_map_struct&);
    void init_headers();
    void print_header(std::string);
    std::unique_ptr<AbstractOntology> spawn_object();

};


#endif //ENTAP_ONTOLOGY_H
