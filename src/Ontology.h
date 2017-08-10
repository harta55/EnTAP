//
// Created by harta on 5/22/17.
//

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

class QuerySequence;

class Ontology {
    typedef std::map<std::string,std::vector<std::string>> go_struct;
    typedef std::map<std::string, QuerySequence> query_map_struct;

public:

    void execute(query_map_struct&,std::string,std::string);
    Ontology(int,std::string,std::string,std::string,std::string,
             boost::program_options::variables_map &, std::string, GraphingManager*);

private:

    struct interpro_struct {
        double _eval;
        std::map<std::string,std::string> _results;
        std::map<std::string,std::vector<std::string>> _go_map;
    };

    const std::string ONTOLOGY_OUT_PATH     = "ontology/";
    const std::string PROCESSED_OUT_DIR     = "processed/";
    const std::string FIGURE_DIR            = "figures/";
    const std::string OUT_UNANNOTATED_NUCL  = "unannotated_sequences.fnn";
    const std::string OUT_UNANNOTATED_PROT  = "unannotated_sequences.faa";
    const std::string OUT_ANNOTATED_NUCL    = "annotated_sequences.fnn";
    const std::string OUT_ANNOTATED_PROT    = "annotated_sequences.faa";
    const std::string GO_MOLECULAR_FLAG     = "molecular_function";
    const std::string GO_BIOLOGICAL_FLAG    = "biological_process";
    const std::string GO_CELLULAR_FLAG      = "cellular_component";
    const std::string GO_OVERALL_FLAG       = "overall";

    const std::string GRAPH_EGG_TAX_BAR_TITLE = "Top_10_Tax_Levels";
    const std::string GRAPH_EGG_TAX_BAR_PNG   = "eggnog_tax_scope.png";
    const std::string GRAPH_EGG_TAX_BAR_TXT   = "eggnog_tax_scope.txt";
    const std::string GRAPH_GO_END_TXT        = "_go_bar_graph.txt";
    const std::string GRAPH_GO_END_PNG        = "_go_bar_graph.png";
    const std::string GRAPH_GO_BAR_BIO_TITLE  = "Top_10_GO_Biological_Terms";
    const std::string GRAPH_GO_BAR_CELL_TITLE = "Top_10_GO_Cellular_Terms";
    const std::string GRAPH_GO_BAR_MOLE_TITLE = "Top_10_GO_Molecular_Terms";
    const std::string GRAPH_GO_BAR_ALL_TITLE  = "Top_10_GO_Terms";
    const unsigned char GRAPH_ONTOLOGY_FLAG = 4;
    const unsigned char GRAPH_TOP_BAR_FLAG= 1;  // used for tax levels and go term tops
    static constexpr short INTERPRO_COL_NUM = 15;
    static constexpr short EGGNOG_COL_NUM   = 12;

    std::vector<short>              _go_levels;
    int                             _threads;
    short                           _software_flag;
    bool                            _is_overwrite;
    std::string                     _entap_exe;
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

    void parse_results_eggnog(query_map_struct&,std::pair<std::string,std::string>&);
    void run_eggnog(query_map_struct&);
    void run_interpro(query_map_struct&,std::vector<std::string>&);
    void parse_results_interpro(query_map_struct&, std::pair<std::string,std::string>&);
    void print_eggnog(query_map_struct&);
    void print_interpro(query_map_struct&);
    go_struct parse_go_list(std::string, std::map<std::string,struct_go_term> &
            ,char);
    std::string eggnog_format(std::string);
    void init_headers();
    void print_header(std::string);
    bool verify_files(std::string,std::string);
    void interpro_format_fix(std::string&);
    std::map<std::string,struct_go_term> read_go_map ();
};


#endif //ENTAP_ONTOLOGY_H
