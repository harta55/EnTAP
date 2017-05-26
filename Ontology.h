//
// Created by harta on 5/22/17.
//

#ifndef ENTAP_ONTOLOGY_H
#define ENTAP_ONTOLOGY_H

#include <iostream>
#include <map>
#include "QuerySequence.h"
#include "DatabaseHelper.h"
class QuerySequence;

class Ontology {
    typedef std::map<std::string, QuerySequence> query_map_struct;
    typedef std::map<std::string,std::vector<std::string>> go_struct;
public:

    void execute(short, query_map_struct&,std::string,std::string);
    Ontology(int, bool,std::string,std::string,std::string,std::string);

private:
    const std::string _GO_BIOLOGICAL_FLAG = "biological_process";
    const std::string _GO_CELLULAR_FLAG = "cellular_component";
    const std::string _GO_MOLECULAR_FLAG = "molecular_function";

    int _threads;
    std::string _entap_exe, _eggnog_exe, _outpath, _new_input, _input_no_hits;
    bool _is_overwrite;
    void parse_results_eggnog(query_map_struct&,std::pair<std::string,std::string>&);
    void run_eggnog(query_map_struct&);
    void print_eggnog(query_map_struct&);
    go_struct parse_go_list(std::string, DatabaseHelper&);
};


#endif //ENTAP_ONTOLOGY_H
