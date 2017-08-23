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

//*********************** Defines ******************************

#define PATHS(x,y)      (boostFS::path(x) / boostFS::path(y)).string()
#define NCBI_UNIPROT    0       // Compiler flag for future feature

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
extern std::string RSEM_EXE_DIR;
extern std::string GENEMARK_EXE;
extern std::string DIAMOND_EXE;
extern std::string EGG_EMAPPER_EXE;
extern std::string EGG_SQL_DB_PATH;
extern std::string EGG_DOWNLOAD_EXE;
extern std::string INTERPRO_EXE;
extern std::string TAX_DB_PATH;
extern std::string TAX_DOWNLOAD_EXE;
extern std::string GO_DB_PATH;
extern std::string GRAPHING_EXE;


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


namespace ENTAP_CONFIG {

    extern const std::string ENTAP_VERSION ;
    extern const std::string DEBUG_FILENAME;
    extern const std::string LOG_FILENAME  ;

    //------------------USER INPUTS-----------------------//
    extern const std::string INPUT_FLAG_CONFIG       ;
    extern const std::string INPUT_FLAG_ALIGN        ;
    extern const std::string INPUT_FLAG_RUNPROTEIN   ;
    extern const std::string INPUT_FLAG_RUNNUCLEOTIDE;
    extern const std::string INPUT_FLAG_OVERWRITE    ;
    extern const std::string INPUT_FLAG_NCBI_1       ;
    extern const std::string INPUT_FLAG_NCBI_2       ;
    extern const std::string INPUT_FLAG_UNIPROT      ;
    extern const std::string INPUT_FLAG_INTERPRO     ;
    extern const std::string INPUT_FLAG_ONTOLOGY     ;
    extern const std::string INPUT_FLAG_SPECIES      ;
    extern const std::string INPUT_FLAG_QCOVERAGE    ;
    extern const std::string INPUT_FLAG_TCOVERAGE    ;
    extern const std::string INPUT_FLAG_COMPLETE     ;
    extern const std::string INPUT_FLAG_GO_LEVELS    ;
    extern const std::string INPUT_FLAG_EXE_PATH     ;
    extern const std::string INPUT_FLAG_FPKM         ;
    extern const std::string INPUT_FLAG_DATA_OUT     ;
    extern const std::string INPUT_FLAG_CONTAM       ;
    extern const std::string INPUT_FLAG_E_VAL        ;
    extern const std::string INPUT_FLAG_HELP         ;
    extern const std::string INPUT_FLAG_VERSION      ;
    extern const std::string INPUT_FLAG_TRANSCRIPTOME;
    extern const std::string INPUT_FLAG_DATABASE     ;

    extern const std::string INPUT_UNIPROT_SWISS    ;
    extern const std::string INPUT_UNIPROT_UR100    ;
    extern const std::string INPUT_UNIPROT_UR90     ;
    extern const std::string INPUT_UNIPROT_TREMBL   ;
    extern const std::string INPUT_UNIPROT_NULL     ;
    extern const std::string INPUT_UNIPROT_DEFAULT  ;

    extern const std::string UNIPROT_BASE_PATH ;
    extern const std::string UNIPROT_INDEX_PATH;

    extern const std::string NCBI_NONREDUNDANT ;
    extern const std::string NCBI_BASE_PATH ;
    extern const std::string NCBI_REFSEQ_COMP ;
    extern const std::string NCBI_REFSEQ_PLANT ;
    extern const std::string NCBI_NULL;
    extern const std::string NCBI_DEFAULT;

    extern const std::string GO_DB_PATH_DEF    ;
    extern const std::string GO_TERM_FILE      ;
    extern const std::string TAX_DATABASE_PATH ;
    extern const std::string TAX_BIN_PATH      ;
    extern const std::string BIN_PATH          ;
    extern const std::string DATABASE_DIR      ;
    extern const std::string NCBI_INDEX_PATH   ;
}

//**************************************************************


//******************* Global Constants *************************
namespace ENTAP_ERR {
    const unsigned short E_INPUT_PARSE = 10;
    const unsigned short E_SUCCESS = 11;
    const unsigned short E_CONFIG_PARSE = 12;
    const unsigned short E_CONFIG_CREATE = 13;
    const unsigned short E_CONFIG_CREATE_SUCCESS = 14;
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

namespace ENTAP_STATS {
    const std::string SOFTWARE_BREAK = "----------------------------------------------\n";
}

//***********************************************




#endif //ENTAPGLOBALS_H
