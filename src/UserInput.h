/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/

#ifndef ENTAP_USERINPUT_H
#define ENTAP_USERINPUT_H

//*********************** Includes *****************************
#include <boost/program_options/variables_map.hpp>
#include "EntapGlobals.h"

//**************************************************************


//*********************** Typedefs *****************************
typedef std::vector<std::string> databases_t;

//******************** Prototype Functions *********************
boost::program_options::variables_map parse_arguments_boost(int, const char**);
bool verify_user_input(boost::program_options::variables_map&);
void print_user_input(boost::program_options::variables_map &map, std::string&, std::string&);
bool check_key(std::string&);
std::unordered_map<std::string,std::string> parse_config(std::string&,std::string&);
void generate_config(std::string&);
void verify_databases(boost::program_options::variables_map&);
void verify_species (boost::program_options::variables_map&);
void init_exe_paths(std::unordered_map<std::string, std::string> &, std::string);
std::string get_exe_path(boostPO::variables_map&);

//**************************************************************


//*********************** Constants ****************************
const float DEFAULT_QCOVERAGE              = 50.0;
const float DEFAULT_TCOVERAGE              = 50.0;
const float E_VALUE                        = 1e-5;
const float RSEM_FPKM_DEFAULT              = 0.5;
const unsigned char MAX_DATABASE_SIZE      = 5;
const std::string INTERPRO_DEFAULT         = "pfam";
const std::string OUTFILE_DEFAULT          = "outfiles";


//-------------------Config File----------------------//
const std::string CONFIG_FILE              = "entap_config.txt";
const std::string KEY_UNIPROT_SWISS        = "uniprot_swiss_path";
const std::string KEY_UNIPROT_UR90         = "uniprot_ur90_path";
const std::string KEY_UNIPROT_UR100        = "uniprot_ur100_path";
const std::string KEY_UNIPROT_TREMBL       = "uniprot_trembl_path";
const std::string KEY_NCBI_NR              = "ncbi_nr_path";
const std::string KEY_NCBI_REFSEQ_COMPLETE = "ncbi_refseq_complete_path";
const std::string KEY_NCBI_REFSEQ_SEPARATE = "ncbi_refseq_separate_path";
const std::string KEY_DIAMOND_EXE          = "diamond_exe_path";
const std::string KEY_RSEM_EXE             = "rsem_exe_path";
const std::string KEY_GENEMARK_EXE         = "genemarkst_exe_path";
const std::string KEY_EGGNOG_EXE           = "eggnog_exe_path";
const std::string KEY_EGGNOG_DOWN          = "eggnog_download_exe";
const std::string KEY_INTERPRO_EXE         = "interpro_exe_path";
const std::string KEY_EGGNOG_DB            = "eggnog_database";
const std::string KEY_TAX_DB               = "entap_tax_database";
const std::string KEY_GO_DB                = "entap_go_database";
const std::string KEY_TAX_DOWNLOAD_EXE     = "entap_tax_download_script";
const std::string KEY_GRAPH_SCRIPT         = "entap_graphing_script";

// Avoid cluttering global namespace / conflicts for config paths
namespace Defaults {
    const std::string RSEM_DEFAULT_EXE         = "/libs/RSEM-1.3.0/";
    const std::string GENEMARK_DEFAULT_EXE     = "/libs/gmst_linux_64/gmst.pl";
    const std::string DIAMOND_DEFAULT_EXE      = "/libs/diamond-0.8.31/bin/diamond";
    const std::string EGG_EMAPPER_DEFAULT      = "/libs/eggnog-mapper/emapper.py";
    const std::string EGG_DOWNLOAD_DEFAULT     = "/libs/eggnog-mapper/download_eggnog_data.py";
    const std::string EGG_SQL_DB_DEFAULT       = "/libs/eggnog-mapper/data/eggnog.db";
    const std::string INTERPRO_DEF_EXE         = "/libs/interproscan-5.22-61.0/interproscan.sh";
    const std::string TAX_DOWNLOAD_DEF         = "/src/download_tax.pl";
    const std::string GRAPH_SCRIPT_DEF         = "/src/entap_graphing.py";

}

//**************************************************************






#endif //ENTAP_USERINPUT_H
