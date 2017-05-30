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

    void execute(short, query_map_struct&,std::string,std::string,std::vector<std::string>&);
    Ontology(int, bool,std::string,std::string,std::string,std::string);

private:
    std::string ONTOLOGY_OUT_PATH;
    std::vector<std::string> _HEADERS;
    int _threads;
    short _software_flag;
    std::string _entap_exe, _ontology_exe, _outpath, _new_input, _input_no_hits;
    bool _is_overwrite;
    void parse_results_eggnog(query_map_struct&,std::pair<std::string,std::string>&);
    void run_eggnog(query_map_struct&);
    void run_interpro(query_map_struct&,std::vector<std::string>&);
    void parse_results_interpro(query_map_struct&, std::pair<std::string,std::string>&);
    void print_eggnog(query_map_struct&);
    void print_interpro(query_map_struct&);
    go_struct parse_go_list(std::string, DatabaseHelper&, char);
    std::string eggnog_format(std::string);
    void init_headers();
    void print_header(std::string);
};


#endif //ENTAP_ONTOLOGY_H
