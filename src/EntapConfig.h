//
// Created by harta55 on 2/1/17.
//

#ifndef ENTAPCONFIG_H
#define ENTAPCONFIG_H

#include <string>
#include <boost/serialization/unordered_map.hpp>
#include <boost/program_options/variables_map.hpp>


namespace entapConfig {
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
    const std::string ENTAP_CONFIG_DIR = "/entap_config";
    const std::string TAX_DATABASE_PATH = "/databases/ncbi_tax.entp";

    //******************Prototype Functions******************

    void init_entap(boost::program_options::variables_map, std::string);
    void init_taxonomic(std::string&);
    void init_uniprot(std::vector<std::string>&, std::string);
    void init_ncbi(std::vector<std::string>&, std::string);
    void init_diamond_index(std::string,std::string,int);
    std::string download_file(std::string, std::string&,std::string&);
    std::string download_file(const std::string &,std::string&);
    void decompress_file(std::string,std::string,short);
    int update_database(std::string);
    void init_go_db(std::string&,std::string);
    void init_eggnog(std::string);


}

#endif //ENTAP_INITHANDLER_H
