//
// Created by harta on 5/22/17.
//

#ifndef ENTAP_ONTOLOGY_H
#define ENTAP_ONTOLOGY_H

#include <iostream>
#include <map>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"
#include "EntapInit.h"

class QuerySequence;

class Ontology {
    typedef std::map<std::string, QuerySequence> query_map_struct;
    typedef std::map<std::string,std::vector<std::string>> go_struct;
public:

    void execute(query_map_struct&,std::string,std::string);
    Ontology(int,std::string,std::string,std::string,std::string,
             boost::program_options::variables_map &);

private:
    const std::string ONTOLOGY_OUT_PATH = "ontology/";
    const std::string PROCESSED_OUT_DIR = "ontology/processed/";
    const std::string OUT_UNANNOTATED_NUCL = "unannotated_sequences.fnn";
    const std::string OUT_UNANNOTATED_PROT = "unannotated_sequences.faa";
    const std::string OUT_ANNOTATED_NUCL = "annotated_sequences.fnn";
    const std::string OUT_ANNOTATED_PROT = "annotated_sequences.faa";
    std::vector<std::string> _HEADERS, _interpro_databases;
    std::vector<short> _go_levels;
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
    go_struct parse_go_list(std::string, std::map<std::string,entapInit::struct_go_term> &
            ,char);
    std::string eggnog_format(std::string);
    void init_headers();
    void print_header(std::string);
    bool verify_files(std::string,std::string);
    void interpro_format_fix(std::string&);
    std::map<std::string,entapInit::struct_go_term> read_go_map ();
};


#endif //ENTAP_ONTOLOGY_H
