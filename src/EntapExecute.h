//
// Created by harta on 3/4/17.
//

#ifndef ENTAP_ENTAPEXECUTE_H
#define ENTAP_ENTAPEXECUTE_H

#include "Ontology.h"

#include "QuerySequence.h"
#include <list>
#include <boost/program_options/variables_map.hpp>
#include <map>
#include <unordered_map>
#include <queue>

namespace entapExecute {

    enum ExecuteStates {
        INIT,
        FRAME_SELECTION,
        RSEM,
        FILTER,
        DIAMOND_RUN,
        DIAMOND_PARSE,
        GENE_ONTOLOGY,
        EXIT
    };

    std::list<std::string> verify_databases(std::vector<std::string>, std::vector<std::string>,
                                            std::vector<std::string>, std::string &,
                                            std::unordered_map<std::string, std::string> &);

    void execute_main(boost::program_options::variables_map &, std::string,
                      std::unordered_map<std::string, std::string> &);
    std::string filter_transcriptome(std::string &, std::string &, float, std::string, bool);
    void print_filtered_map(std::map<std::string, QuerySequence> &, std::string &);
    void verify_state(std::queue<char> &, bool &);
    bool is_file_empty(std::string);
    bool valid_state(enum ExecuteStates);
    std::string init_exe_paths(std::unordered_map<std::string, std::string> &, std::string &);
    void print_statistics(std::string &, std::string &);
    std::map<std::string, QuerySequence> init_sequence_map(std::string&);
    std::pair<unsigned long, unsigned long> calculate_N_vals
            (std::vector<unsigned long> &, unsigned long);
}


#endif //ENTAP_ENTAPEXECUTE_H
