/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
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

class Ontology {
    typedef std::map<std::string, QuerySequence> query_map_struct;


public:

    void execute(query_map_struct&,std::string,std::string);
    Ontology(int,std::string,std::string, boost::program_options::variables_map &,GraphingManager*);

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
    GraphingManager  *_graphingManager;

    void print_eggnog(query_map_struct&);
    void print_interpro(query_map_struct&);
    void init_headers();
    void print_header(std::string);
    std::unique_ptr<AbstractOntology> spawn_object();

};


#endif //ENTAP_ONTOLOGY_H
