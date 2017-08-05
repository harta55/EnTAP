//
// Created by harta on 5/10/17.
//

#ifndef ENTAP_SIMILARITYSEARCH_H
#define ENTAP_SIMILARITYSEARCH_H

//*********************** Includes *****************************
#include <iostream>
#include <list>
#include <unordered_map>
#include <map>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"
#include "GraphingManager.h"

//**************************************************************


class SimilaritySearch {


public:

    //******************** Public Prototype Functions *********************
    std::vector<std::string> execute(std::string, bool);
    SimilaritySearch(std::vector<std::string>&, std::string, int, std::string, std::string,
                     std::string,boost::program_options::variables_map &, GraphingManager *);
    std::pair<std::string,std::string> parse_files(std::string,
                                                   std::map<std::string, QuerySequence>&);
    //**************************************************************


private:


    const std::string _NCBI_REGEX                                = "\\[([^]]+)\\](?!.+\\[.+\\])";
    const std::string _UNIPROT_REGEX                             = "OS=(.+?)\\s\\S\\S=";
    const std::string SIM_SEARCH_DATABASE_BEST_TSV               = "best_hits.tsv";
    const std::string SIM_SEARCH_DATABASE_BEST_TSV_NO_CONTAM     = "best_hits_no_contam.tsv";
    const std::string SIM_SEARCH_DATABASE_BEST_FA_NUCL           = "best_hits.fnn";
    const std::string SIM_SEARCH_DATABASE_BEST_FA_PROT           = "best_hits.faa";
    const std::string SIM_SEARCH_DATABASE_BEST_FA_NUCL_NO_CONTAM = "best_hits_no_contam.fnn";
    const std::string SIM_SEARCH_DATABASE_BEST_FA_PROT_NO_CONTAM = "best_hits_no_contam.faa";

    const std::string SIM_SEARCH_DATABASE_CONTAM_TSV             = "best_hits_contam.tsv";
    const std::string SIM_SEARCH_DATABASE_CONTAM_FA_NUCL         = "best_hits_contam.fnn";
    const std::string SIM_SEARCH_DATABASE_CONTAM_FA_PROT         = "best_hits_contam.faa";
    const std::string SIM_SEARCH_DATABASE_NO_HITS_NUCL           = "no_hits.fnn";
    const std::string SIM_SEARCH_DATABASE_NO_HITS_PROT           = "no_hits.faa";
    const std::string SIM_SEARCH_DATABASE_UNSELECTED             = "unselected.tsv";
    const std::string SIM_SEARCH_DIR                             = "similarity_search/";
    const std::string PROCESSED_DIR                              = "processed/";
    const std::string RESULTS_DIR                                = "overall_results/";
    const std::string FIGURE_DIR                                 = "figures/";

    const unsigned char GRAPH_SOFTWARE_FLAG                      = 3;
    const unsigned char GRAPH_BAR_FLAG                           = 1;
    std::string GRAPH_SPECIES_BAR_TXT                            = "_species_bar.txt";
    std::string GRAPH_SPECIES_BAR_PNG                            = "_species_bar.png";
    std::string GRAPH_SPECIES_TITLE                              = "_Top_10_Species_Distribution";
    std::string GRAPH_CONTAM_BAR_TXT                             = "_contam_bar.txt";
    std::string GRAPH_CONTAM_BAR_PNG                             = "_contam_bar.png";
    std::string GRAPH_CONTAM_TITLE                               = "_Top_10_Contaminant_Distribution";

    static constexpr int DMND_COL_NUMBER = 14;

    const std::list<std::string> INFORMATIVENESS {
            "conserved",
            "predicted",
            "unnamed",
            "hypothetical",
            "putative",
            "unidentified",
            "uncharacterized",
            "unknown",
            "uncultured",
            "uninformative"
    };

    std::vector<std::string>         _database_paths;
    std::vector<std::string>         _sim_search_paths;
    std::string                      _diamond_exe;
    std::string                      _outpath;
    std::string                      _input_path;
    std::string                      _sim_search_dir;
    std::string                      _processed_path;
    std::string                      _figure_path;
    std::string                      _results_path;
    std::string                      _entap_exe;
    std::string                      _input_lineage;
    std::string                      _input_species;
    std::string                      _blast_type;
    int                              _threads;
    bool                             _overwrite;
    bool                             _blastp;
    float                           _e_val;
    float                           _qcoverage;
    float                           _tcoverage;
    short                            _software_flag;
    std::vector<std::string>         _contaminants;
    GraphingManager                  *_pGraphingManager;
    std::unordered_map<std::string,std::string> _file_to_database;

    std::vector<std::string> diamond();
    void diamond_blast(std::string, std::string, std::string,std::string&,int&, std::string&);
    std::vector<std::string> verify_diamond_files(std::string&, std::string);
    std::pair<std::string,std::string> diamond_parse(std::vector<std::string>&,
                                                     std::map<std::string, QuerySequence>&);
    std::unordered_map<std::string, std::string> read_tax_map();
    std::pair<bool,std::string>  is_contaminant(std::string, std::unordered_map<std::string,
            std::string> &,std::vector<std::string>&);
    std::string get_species(std::string &);
    bool is_informative(std::string);
    std::pair<std::string,std::string> process_best_diamond_hit(std::list<std::map<std::string,
            QuerySequence>>&,std::map<std::string, QuerySequence>&);
    void print_header(std::string);
    std::string get_lineage(std::string, std::unordered_map<std::string, std::string>&);
    std::pair<std::string,std::string> calculate_best_stats (std::map<std::string, QuerySequence>&,
                                 std::map<std::string, QuerySequence>&,
                                 std::stringstream &, std::string&,bool);
};


#endif //ENTAP_SIMILARITYSEARCH_H
