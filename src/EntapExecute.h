//
// Created by harta on 3/4/17.
//

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
    const std::string RSEM_EXE_PATH        = "/libs/RSEM-1.3.0/";
    const std::string GENEMARK_EXE_PATH    = "/libs/gmst_linux_64/gmst.pl";
    const std::string DIAMOND_PATH_EXE     = "/libs/diamond-0.8.31/bin/diamond";
    const std::string EGGNOG_EMAPPER_EXE   = "/libs/eggnog-mapper/emapper.py";
    const std::string EGGNOG_DOWNLOAD_EXE  = "/libs/eggnog-mapper/download_eggnog_data.py";
    const std::string INTERPRO_EXE         = "/libs/interproscan-5.22-61.0/interproscan.sh";
    const std::string ENTAP_OUTPUT         = "entap_out/";
    const std::string GRAPH_FILEPATH       = "/src/entap_graphing.py";

    //**************************************************************


    // *******************Prototype Functions******************
    std::vector<std::string> verify_databases(std::vector<std::string>, std::vector<std::string>,
                                            std::vector<std::string>, std::string &,
                                            std::unordered_map<std::string, std::string> &);

    void execute_main(boost::program_options::variables_map &, std::string,
                      std::unordered_map<std::string, std::string> &);
    std::string filter_transcriptome(std::string &);
    void verify_state(std::queue<char> &, bool &);
    bool valid_state(enum ExecuteStates);
    std::pair<std::string,std::string> init_exe_paths(std::unordered_map<std::string, std::string> &,
                      std::string);
    std::map<std::string, QuerySequence> init_sequence_map(std::string&,bool);
    std::pair<unsigned long, unsigned long> calculate_N_vals
            (std::vector<unsigned long> &, unsigned long);
    void final_statistics(std::map<std::string, QuerySequence>&);
    //**************************************************************

}


#endif //ENTAP_ENTAPEXECUTE_H
