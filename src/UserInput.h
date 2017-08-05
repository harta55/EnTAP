//
// Created by harta on 8/4/17.
//

#ifndef ENTAP_USERINPUT_H
#define ENTAP_USERINPUT_H

//*********************** Includes *****************************
#include <boost/program_options/variables_map.hpp>

//**************************************************************


//*********************** Typedefs *****************************
typedef std::vector<std::string> databases_t;


//******************** Prototype Functions *********************
boost::program_options::variables_map parse_arguments_boost(int, const char**);
void verify_user_input(boost::program_options::variables_map&);
void print_user_input(boost::program_options::variables_map &map, std::string&, std::string&);
bool check_key(std::string&);
std::unordered_map<std::string,std::string> parse_config(std::string&);
void generate_config(std::string&);
void verify_databases(boost::program_options::variables_map&);
void verify_species (boost::program_options::variables_map&);
std::pair<bool, boost::program_options::variables_map> entap_user_parse(int argc, const char** argv);
//**************************************************************


//*********************** Constants ****************************
const float DEFAULT_QCOVERAGE              = 50.0;
const float DEFAULT_TCOVERAGE              = 50.0;
const float E_VALUE                        = 1e-5;
const float RSEM_FPKM_DEFAULT              = 0.5;
const unsigned char MAX_DATABASE_SIZE      = 5;
const std::string OUTFILE_DEFAULT          = "outfiles";
const std::string INTERPRO_DEFAULT         = "pfam";


//**************************************************************






#endif //ENTAP_USERINPUT_H
