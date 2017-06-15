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
    std::string download_file(const std::string &,std::string&);
    void decompress_file(std::string,std::string,short);
    int update_database(std::string);
    int get_supported_threads(boost::program_options::variables_map&);
    void init_go_db(std::string&,std::string);
    static int callback(void *data, int argc, char **argv, char **azColName);
    std::string generate_command(std::unordered_map<std::string,std::string>&,
                                 std::string);
}

#endif //ENTAP_INITHANDLER_H
