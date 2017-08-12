//
// Created by harta55 on 2/1/17.
//
// TODO temporary, change

#ifndef ENTAPGLOBALS_H
#define ENTAPGLOBALS_H

//*********************** Includes *****************************
#include <vector>
#include <list>
#include <boost/serialization/access.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <unordered_map>

namespace boostFS = boost::filesystem;
namespace boostPO = boost::program_options;

//**************************************************************

//***************** Global Prototype Functions *****************
void print_debug(std::string);
void print_statistics(std::string &msg);
bool file_exists (std::string);
int execute_cmd(std::string,std::string);
int execute_cmd(std::string);
std::string generate_command(std::unordered_map<std::string,std::string>&,
                             std::string);
int get_supported_threads(boost::program_options::variables_map&);

//**************************************************************


//**************** Global Structures/Typedefs ******************
typedef struct {
    std::string     text_file_path;
    std::string     graph_title;
    std::string     fig_out_path;
    unsigned char   software_flag;
    unsigned char   graph_type;
} GraphingStruct;

struct  struct_go_term {
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

typedef std::pair<std::string,int> count_pair;
struct compair {
    bool operator ()(count_pair const& one, count_pair const& two) const {
        return one.second > two.second;
    }
};


//*********************** Externs *****************************

extern std::string DEBUG_FILE_PATH;
extern std::string LOG_FILE_PATH;


namespace ENTAP_EXECUTE {
    //------------------------Ontology-------------------------//
    extern const std::string GO_BIOLOGICAL_FLAG ;
    extern const std::string GO_CELLULAR_FLAG;
    extern const std::string GO_MOLECULAR_FLAG;
    const short EGGNOG_INT_FLAG = 0;
    const short INTERPRO_INT_FLAG = 1;


    //------------------------Headers-------------------------//
    extern const std::string HEADER_QUERY;
    extern const std::string HEADER_SUBJECT;
    extern const std::string HEADER_PERCENT;
    extern const std::string HEADER_ALIGN_LEN;
    extern const std::string HEADER_MISMATCH;
    extern const std::string HEADER_GAP_OPEN;
    extern const std::string HEADER_QUERY_S;
    extern const std::string HEADER_QUERY_E;
    extern const std::string HEADER_SUBJ_S;
    extern const std::string HEADER_SUBJ_E;
    extern const std::string HEADER_E_VAL;
    extern const std::string HEADER_COVERAGE;
    extern const std::string HEADER_TITLE;
    extern const std::string HEADER_SPECIES;
    extern const std::string HEADER_DATABASE;
    extern const std::string HEADER_FRAME;
    extern const std::string HEADER_CONTAM;

    extern const std::string HEADER_SEED_ORTH;
    extern const std::string HEADER_SEED_EVAL;
    extern const std::string HEADER_SEED_SCORE;
    extern const std::string HEADER_PRED_GENE;
    extern const std::string HEADER_TAX_SCOPE;
    extern const std::string HEADER_EGG_OGS;
    extern const std::string HEADER_EGG_KEGG;
    extern const std::string HEADER_EGG_GO_BIO ;
    extern const std::string HEADER_EGG_GO_CELL;
    extern const std::string HEADER_EGG_GO_MOLE;
    extern const std::string HEADER_EGG_DESC;
    extern const std::string HEADER_EGG_LEVEL;
    extern const std::string HEADER_EGG_PROTEIN;
}

//**************************************************************


//******************* Global Constants *************************
namespace ENTAP_ERR {
    const unsigned short E_INPUT_PARSE = 10;
    const unsigned short E_SUCCESS = 11;
    const unsigned short E_CONFIG_PARSE = 12;
    const unsigned short E_CONFIG_CREATE = 13;
    const unsigned short E_INIT_TAX_DOWN = 20;
    const unsigned short E_INIT_TAX_INDEX = 21;
    const unsigned short E_INIT_TAX_SERIAL = 22;
    const unsigned short E_INIT_INDX_DATA_NOT_FOUND = 30;
    const unsigned short E_INIT_INDX_DATABASE = 31;
    const unsigned short E_INIT_DOWNLOAD = 23;
    const unsigned short E_INIT_EGGNOG   = 40;

    const unsigned short E_INIT_TAX_READ = 55;
    const unsigned short E_INIT_GO_SETUP = 60;

    const unsigned short E_RUN_EXECUTION_PATHS         = 105;
    const unsigned short E_RUN_VERIFY_DATABASES        = 106;
    const unsigned short E_RUN_GENEMARK                = 100;
    const unsigned short E_RUN_GENEMARK_PARSE          = 101;
    const unsigned short E_RUN_GENEMARK_STATS          = 102;
    const unsigned short E_RUN_RSEM_VALIDATE           = 110;
    const unsigned short E_RUN_RSEM_CONVERT            = 111;
    const unsigned short E_RUN_RSEM_EXPRESSION         = 112;
    const unsigned short E_RUN_FILTER                  = 120;
    const unsigned short E_RUN_SIM_SEARCH_FILTER       = 140;
    const unsigned short E_RUN_ANNOTATION              = 150;
    const unsigned short E_RUN_EGGNOG                  = 160;
    const unsigned short E_PARSE_EGGNOG                = 170;
}

namespace ENTAP_CONFIG {

    const std::string ENTAP_VERSION  = "0.7.1.1";
    const std::string DEBUG_FILENAME = "debug.txt";
    const std::string LOG_FILENAME   = "log_file.txt";

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
    const std::string KEG_EGGNOG_DMND_DB       = "eggnog_dmnd_database";

    //------------------USER INPUTS-----------------------//
    const std::string INPUT_FLAG_CONFIG        = "config";
    const std::string INPUT_FLAG_ALIGN         = "align";
    const std::string INPUT_FLAG_RUNPROTEIN    = "runP";
    const std::string INPUT_FLAG_RUNNUCLEOTIDE = "runN";
    const std::string INPUT_FLAG_OVERWRITE     = "overwrite";
    const std::string INPUT_FLAG_NCBI_1        = "ncbi";
    const std::string INPUT_FLAG_NCBI_2        = "N";
    const std::string INPUT_FLAG_UNIPROT       = "uniprot";
    const std::string INPUT_FLAG_INTERPRO      = "protein";
    const std::string INPUT_FLAG_ONTOLOGY      = "ontology";
    const std::string INPUT_FLAG_SPECIES       = "species";
    const std::string INPUT_FLAG_QCOVERAGE     = "qcoverage";
    const std::string INPUT_FLAG_TCOVERAGE     = "tcoverage";
    const std::string INPUT_FLAG_COMPLETE      = "complete";
    const std::string INPUT_FLAG_GO_LEVELS     = "level";
    const std::string INPUT_FLAG_EXE_PATH      = "exe";
    const std::string INPUT_FLAG_FPKM          = "fpkm";
    const std::string INPUT_FLAG_DATA_OUT      = "database-out";
    const std::string INPUT_FLAG_CONTAM        = "contam";
    const std::string INPUT_FLAG_E_VAL         = "e";
    const std::string INPUT_FLAG_HELP          = "help";
    const std::string INPUT_FLAG_VERSION       = "version";
    const std::string INPUT_FLAG_TRANSCRIPTOME = "input";
    const std::string INPUT_FLAG_DATABASE      = "database";

    const std::string INPUT_UNIPROT_SWISS      = "swiss";
    const std::string INPUT_UNIPROT_UR100      = "ur100";
    const std::string INPUT_UNIPROT_UR90       = "ur90";
    const std::string INPUT_UNIPROT_TREMBL     = "trembl";
    const std::string INPUT_UNIPROT_NULL       = "null";
    const std::string INPUT_UNIPROT_DEFAULT    = INPUT_UNIPROT_SWISS;

    const std::string UNIPROT_BASE_PATH = "/databases/uniprot_";
    const std::string UNIPROT_INDEX_PATH = "/bin/uniprot_";

    const std::string NCBI_NONREDUNDANT = "nr";
    const std::string NCBI_BASE_PATH = "/databases/ncbi_";
    const std::string NCBI_REFSEQ_COMP = "refseq-c";
    const std::string NCBI_REFSEQ_PLANT = "refseq-p";
    const std::string NCBI_NULL = "null";
    const std::string NCBI_DEFAULT = NCBI_REFSEQ_COMP;

    const std::string GO_DB_PATH        = "/bin/go_term.entp";
    const std::string GO_TERM_FILE      = "/databases/term.txt";
    const std::string TAX_SCRIPT_PATH   = "/src/download_tax.pl";
    const std::string TAX_DATABASE_PATH = "/databases/ncbi_tax.entp";
    const std::string TAX_BIN_PATH      = "/bin/ncbi_tax_bin.entp";
    const std::string BIN_PATH          = "bin/";
    const std::string DATABASE_DIR      = "databases/";
    const std::string NCBI_INDEX_PATH   = "/bin/ncbi_";
}


namespace ENTAP_STATS {
    const std::string SOFTWARE_BREAK = "----------------------------------------------\n";
}

//***********************************************




#endif //ENTAPGLOBALS_H
