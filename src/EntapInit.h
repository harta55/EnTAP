//
// Created by harta55 on 2/1/17.
//

#ifndef ENTAP_INITHANDLER_H
#define ENTAP_INITHANDLER_H

#include <string>
#include <boost/serialization/unordered_map.hpp>
#include <boost/program_options/variables_map.hpp>


namespace entapInit {
    struct struct_go_term {
        std::string go_id, level, category, term;
        friend class boost::serialization::access;
        template <typename Archive>
        void serialize(Archive & ar, const unsigned int v) {
            ar&go_id;
            ar&level;
            ar&category;
            ar&term;
        }
    };

    //-----------------------FTP PATHS---------------------------//
    const std::string GO_DATABASE_FTP =
            "http://archive.geneontology.org/latest-full/go_monthly-termdb-tables.tar.gz";
    const std::string UNIPROT_FTP_SWISS = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";
    const std::string UNIPROT_FTP_TREMBL = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz";
    const std::string GO_TERM_FILE = "term.txt";
    const std::string GO_GRAPH_FILE = "graph_path.txt";
    const std::string GO_DATA_NAME = "go_monthly-termdb-tables.tar.gz";
    const std::string GO_DIR = "go_monthly-termdb-tables/";

    /******************Prototype Functions******************/

    bool file_exists (const std::string& name);
    void print_msg(std::string msg);
    void init_entap(boost::program_options::variables_map, std::string,
        std::unordered_map<std::string,std::string>&);
    void init_taxonomic(std::string&);
    void init_uniprot(std::vector<std::string>&, std::string);
    void init_ncbi(std::vector<std::string>&, std::string);
    void init_diamond_index(std::string,std::string,int);
    int execute_cmd(std::string,std::string);
    int execute_cmd(std::string);
    std::string download_file(std::string, std::string&,std::string&);
    std::string download_file(const std::string &,std::string&);
    void decompress_file(std::string,std::string,short);
    int update_database(std::string);
    int get_supported_threads(boost::program_options::variables_map&);
    void init_go_db(std::string&,std::string);
    std::string generate_command(std::unordered_map<std::string,std::string>&,
                                 std::string);
    void init_eggnog(std::string);
}

#endif //ENTAP_INITHANDLER_H
