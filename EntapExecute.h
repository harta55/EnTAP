//
// Created by harta on 3/4/17.
//

#ifndef ENTAP_ENTAPEXECUTE_H
#define ENTAP_ENTAPEXECUTE_H


#include "QuerySequence.h"

namespace entapExecute {
    void execute_main(std::unordered_map<std::string, std::string>);
    void genemarkST();
    void rsem();
    void diamond_run(std::string, std::string, std::string);
    void diamond_parse(std::string[], int);
    void diamond_blast(std::string ,std::string, std::string);
    bool is_contaminant(std::string);
    void print_map(std::unordered_map<std::string, QuerySequence>);
    std::unordered_map<std::string, std::string> read_tax_map();

};


#endif //ENTAP_ENTAPEXECUTE_H
