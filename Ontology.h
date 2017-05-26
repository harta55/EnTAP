//
// Created by harta on 5/22/17.
//

#ifndef ENTAP_ONTOLOGY_H
#define ENTAP_ONTOLOGY_H


#include <iostream>
#include <map>
#include "QuerySequence.h"

class Ontology {
    typedef std::map<std::string, QuerySequence> query_map_struct;

public:
    void execute(short, query_map_struct&,std::string,std::string);
    Ontology(int, bool,std::string,std::string,std::string);

private:
    int _threads;
    std::string _entap_exe, _eggnog_exe, _outpath, _new_input, _input_no_hits;
    bool _is_overwrite;
    void parse_results_eggnog(query_map_struct&,std::pair<std::string,std::string>&);
    void run_eggnog(std::map<std::string,QuerySequence>&);
    void print_eggnog(query_map_struct&);
};


#endif //ENTAP_ONTOLOGY_H
