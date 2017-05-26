//
// Created by harta55 on 2/1/17.
//

#ifndef ENTAP_INITHANDLER_H
#define ENTAP_INITHANDLER_H

#include <string>
#include <boost/serialization/unordered_map.hpp>
#include <boost/program_options/variables_map.hpp>


namespace entapInit {

    struct FtpFile {
        const char *filename;
        FILE *stream;
    };

    bool file_exists (const std::string& name);
    void print_msg(std::string msg);
    void init_entap(boost::program_options::variables_map, std::string,
        std::unordered_map<std::string,std::string>&);
    void init_taxonomic(std::string&);
    void init_uniprot(std::vector<std::string>&, std::string);
    void init_ncbi(std::vector<std::string>&, std::string);
    void init_diamond_index(std::string,std::string,int);
    void verify_state();
    int execute_cmd(std::string,std::string);
    int execute_cmd(std::string);
    std::string download_file(std::string, std::string&,std::string&);
    std::string download_file(std::string &,std::string&);
    void decompress_file(std::string);
    int update_database(std::string);
    int get_supported_threads(boost::program_options::variables_map&);
    void init_go_db(std::string&);
    static int callback(void *data, int argc, char **argv, char **azColName);
}

#endif //ENTAP_INITHANDLER_H
