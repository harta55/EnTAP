//
// Created by harta on 5/10/17.
//

#ifndef ENTAP_SIMILARITYSEARCH_H
#define ENTAP_SIMILARITYSEARCH_H
#include <iostream>
#include <list>
#include <unordered_map>
#include <map>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"


class SimilaritySearch {
public:
    std::list<std::string> execute(std::string, bool);
    SimilaritySearch(std::list<std::string>&, std::string, int, std::string, std::string,
                     std::string,boost::program_options::variables_map &);
    std::pair<std::string,std::string> parse_files(std::string,std::map<std::string, QuerySequence>&);
private:
    std::list<std::string> _database_paths, _sim_search_paths;
    std::string _diamond_exe, _outpath, _input_path, _entap_exe, _input_lineage,_input_species;
    std::string _ncbi_regex = "\\[([^]]+)\\](?!.+\\[.+\\])";
    std::string _uniprot_regex = "OS=(.+?)\\s\\S\\S=";
    int _threads; bool _overwrite,_blastp;
    double _e_val, _coverage;
    short _software_flag;
    std::vector<std::string> _contaminants;
    std::list<std::string> diamond();
    void diamond_blast(std::string, std::string, std::string,std::string&,int&, std::string&);
    std::list<std::string> verify_diamond_files(std::string&,
        std::string&, std::string);
    std::pair<std::string,std::string> diamond_parse(std::vector<std::string>&,
                                                     std::map<std::string, QuerySequence>&);
    std::list<std::string> find_diamond_files();
    std::unordered_map<std::string, std::string> read_tax_map();
    std::pair<bool,std::string>  is_contaminant(std::string, std::unordered_map<std::string,
            std::string> &,std::vector<std::string>&);
    std::string get_species(std::string &);
    bool is_informative(std::string);
    std::pair<std::string,std::string> process_best_diamond_hit(std::list<std::map<std::string,
            QuerySequence>>&,std::map<std::string, QuerySequence>&);
    void print_header(std::string);
    std::string get_lineage(std::string, std::unordered_map<std::string, std::string>&);
    int calculate_score(std::string,bool);
    std::pair<std::string,std::string> calculate_best_stats (std::map<std::string, QuerySequence>&,
                                 std::map<std::string, QuerySequence>&,
                                 std::stringstream &, std::string&,bool);
};


#endif //ENTAP_SIMILARITYSEARCH_H
