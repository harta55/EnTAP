/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/


#ifndef ENTAP_ABSTRACTONTOLOGY_H
#define ENTAP_ABSTRACTONTOLOGY_H

#include <string>
#include "../GraphingManager.h"
#include "../QuerySequence.h"

class QuerySequence;

class AbstractOntology {
public:
    AbstractOntology(std::string &exe,std::string &out, std::string &in_hits,
                     std::string &in_nohits, std::string &proc,
                  std::string &fig, std::string &ont_out, GraphingManager *graphing){
        _exe_path = exe;
        _outpath = out;
        _inpath = in_hits;
        _in_no_hits = in_nohits;
        _processed_path = proc;
        _figure_path = fig;
        _ontology_dir = ont_out;
        pGraphingManager = graphing;

    }
    virtual ~AbstractOntology() = default;
    virtual std::pair<bool, std::string> verify_files()=0;
    virtual void execute(std::map<std::string, QuerySequence>&) = 0;
    virtual void parse(std::map<std::string, QuerySequence>&) = 0;
    virtual void set_data(std::vector<short>&, std::string&, int)=0;

protected:

    const std::string GO_MOLECULAR_FLAG     = "molecular_function";
    const std::string GO_BIOLOGICAL_FLAG    = "biological_process";
    const std::string GO_CELLULAR_FLAG      = "cellular_component";
    const std::string GO_OVERALL_FLAG       = "overall";

    const std::string OUT_UNANNOTATED_NUCL  = "unannotated_sequences.fnn";
    const std::string OUT_UNANNOTATED_PROT  = "unannotated_sequences.faa";
    const std::string OUT_ANNOTATED_NUCL    = "annotated_sequences.fnn";
    const std::string OUT_ANNOTATED_PROT    = "annotated_sequences.faa";

    const std::string GRAPH_GO_END_TXT        = "_go_bar_graph.txt";
    const std::string GRAPH_GO_END_PNG        = "_go_bar_graph.png";
    const std::string GRAPH_GO_BAR_BIO_TITLE  = "Top_10_GO_Biological_Terms";
    const std::string GRAPH_GO_BAR_CELL_TITLE = "Top_10_GO_Cellular_Terms";
    const std::string GRAPH_GO_BAR_MOLE_TITLE = "Top_10_GO_Molecular_Terms";
    const std::string GRAPH_GO_BAR_ALL_TITLE  = "Top_10_GO_Terms";

    const unsigned char GRAPH_ONTOLOGY_FLAG = 4;
    const unsigned char GRAPH_TOP_BAR_FLAG= 1;  // used for tax levels and go term tops

    std::string _in_no_hits;
    std::string _exe_path;
    std::string _outpath;
    std::string _inpath;
    std::string _processed_path;
    std::string _figure_path;
    std::string _ontology_dir;
    GraphingManager *pGraphingManager;

    typedef std::map<std::string,std::vector<std::string>> go_struct;
    typedef std::map<std::string, QuerySequence> query_map_struct;

    std::map<std::string,std::vector<std::string>>
    parse_go_list(std::string list, std::map<std::string,struct_go_term> &GO_DATABASE,char delim);
    std::map<std::string,struct_go_term> read_go_map ();

};


#endif //ENTAP_ABSTRACTONTOLOGY_H
