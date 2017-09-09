/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/

#ifndef ENTAP_ENTAPEXECUTE_H
#define ENTAP_ENTAPEXECUTE_H

//*********************** Includes *****************************
#include "Ontology.h"

#include "QuerySequence.h"
#include <list>
#include <boost/program_options/variables_map.hpp>
#include <map>
#include <unordered_map>
#include <queue>

//**************************************************************


namespace entapExecute {

    enum ExecuteStates {
        INIT,
        RSEM,
        FRAME_SELECTION,
        FILTER,
        DIAMOND_RUN,
        DIAMOND_PARSE,
        GENE_ONTOLOGY,
        EXIT
    };


    // ********************** Global Constants *********************
    const std::string OUT_UNANNOTATED_NUCL = "final_unannotated.fnn";
    const std::string OUT_UNANNOTATED_PROT = "final_unannotated.faa";
    const std::string OUT_ANNOTATED_NUCL   = "final_annotated.fnn";
    const std::string OUT_ANNOTATED_PROT   = "final_annotated.faa";
    const std::string ENTAP_OUTPUT         = "entap_out/";

    //**************************************************************


    // *******************Prototype Functions******************
    std::vector<std::string> verify_databases(std::vector<std::string>, std::vector<std::string>,
                                            std::vector<std::string>, std::string);
    void execute_main(boost::program_options::variables_map &);
    std::string filter_transcriptome(std::string &);
    void verify_state(std::queue<char> &, bool &);
    bool valid_state(enum ExecuteStates);
    std::map<std::string, QuerySequence> init_sequence_map(std::string&,bool,bool);
    std::pair<unsigned long, unsigned long> calculate_N_vals
            (std::vector<unsigned long> &, unsigned long);
    void final_statistics(std::map<std::string, QuerySequence>&);
    void flag_transcripts(ExecuteStates, std::map<std::string, QuerySequence>&);
    //**************************************************************

}


#endif //ENTAP_ENTAPEXECUTE_H
