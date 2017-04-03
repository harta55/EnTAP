//
// Created by harta55 on 2/1/17.
//

#ifndef ENTAP_INITHANDLER_H
#define ENTAP_INITHANDLER_H

#include <string>
#include <boost/serialization/unordered_map.hpp>


namespace entapInit {

    struct FtpFile {
        const char *filename;
        FILE *stream;
    };

    void download_file(std::string, std::string);
    bool file_exists (const std::string& name);
    void print_input(std::unordered_map<std::string, std::string>);
    void print_msg(std::string msg);
    void init_entap(std::unordered_map<std::string, std::string>, std::string);
    void init_taxonomic();
    void init_uniprot(std::string);
    void init_ncbi(std::string);
    void init_database_parse(std::string);
    void init_diamond_index(std::string, std::string, std::string);
    void verify_state();
    int execute_cmd(std::string,std::string);
    int execute_cmd(std::string);
}

#endif //ENTAP_INITHANDLER_H
