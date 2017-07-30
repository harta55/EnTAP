//
// Created by harta55 on 2/1/17.
//

#ifndef ENTAP_ERRORFLAGS_H
#define ENTAP_ERRORFLAGS_H
#include <vector>
#include <list>

// TODO temporary, change

typedef struct {
    std::string     text_file_path;
    std::string     graph_title;
    std::string     fig_out_path;
    unsigned char   software_flag;
    unsigned char   graph_type;
} GraphingStruct;

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

    const std::string ENTAP_VERSION = "0.6.1.1";
    const double DEFAULT_QCOVERAGE = 50.0;
    const double DEFAULT_TCOVERAGE = 50.0;
    const double E_VALUE = 1e-5;
    const std::string DEBUG_FILENAME = "debug.txt";
    const std::string LOG_FILENAME = "log_file.txt";

    //-------------------Config File----------------------//
    const std::string CONFIG_FILE = "entap_config.txt";
    const std::string KEY_UNIPROT_SWISS = "uniprot_swiss_path";
    const std::string KEY_UNIPROT_UR90 = "uniprot_ur90_path";
    const std::string KEY_UNIPROT_UR100 = "uniprot_ur100_path";
    const std::string KEY_UNIPROT_TREMBL = "uniprot_trembl_path";
    const std::string KEY_NCBI_NR = "ncbi_nr_path";
    const std::string KEY_NCBI_REFSEQ_COMPLETE = "ncbi_refseq_complete_path";
    const std::string KEY_NCBI_REFSEQ_SEPARATE = "ncbi_refseq_separate_path";
    const std::string KEY_DIAMOND_EXE = "diamond_exe_path";
    const std::string KEY_RSEM_EXE = "rsem_exe_path";
    const std::string KEY_GENEMARK_EXE = "genemarkst_exe_path";
    const std::string KEY_EGGNOG_EXE = "eggnog_exe_path";
    const std::string KEY_INTERPRO_EXE = "interpro_exe_path";

    //------------------USER INPUTS-----------------------//
    const std::string INPUT_FLAG_ALIGN              = "align";
    const std::string INPUT_FLAG_RUNPROTEIN         = "runP";
    const std::string INPUT_FLAG_RUNNUCLEOTIDE      = "runN";
    const std::string INPUT_FLAG_OVERWRITE          = "overwrite";
    const std::string INPUT_FLAG_NCBI_1             = "ncbi";
    const std::string INPUT_FLAG_NCBI_2             = "N";
    const std::string INPUT_FLAG_UNIPROT            = "uniprot";
    const std::string INPUT_FLAG_INTERPRO           = "protein";
    const std::string INPUT_FLAG_ONTOLOGY           = "ontology";
    const std::string INPUT_FLAG_SPECIES            = "species";
    const std::string INPUT_FLAG_QCOVERAGE          = "qcoverage";
    const std::string INPUT_FLAG_TCOVERAGE          = "tcoverage";
    const std::string INPUT_FLAG_COMPLETE           = "complete";
    const std::string INPUT_FLAG_GO_LEVELS          = "level";
    const std::string INPUT_FLAG_EXE_PATH           = "exe";
    const std::string INPUT_FLAG_FPKM               = "fpkm";

    const std::string INPUT_UNIPROT_SWISS = "swiss";
    const std::string INPUT_UNIPROT_UR100 = "ur100";
    const std::string INPUT_UNIPROT_UR90 = "ur90";
    const std::string INPUT_UNIPROT_TREMBL = "trembl";
    const std::string INPUT_UNIPROT_NULL = "null";
    const std::string INPUT_UNIPROT_DEFAULT = INPUT_UNIPROT_SWISS;

    const std::string UNIPROT_BASE_PATH = "/databases/uniprot_";
    const std::string UNIPROT_INDEX_PATH = "/bin/uniprot_";

    const std::string NCBI_NONREDUNDANT = "nr";
    const std::string NCBI_BASE_PATH = "/databases/ncbi_";
    const std::string NCBI_REFSEQ_COMP = "refseq-c";
    const std::string NCBI_REFSEQ_PLANT = "refseq-p";
    const std::string NCBI_NULL = "null";
    const std::string NCBI_DEFAULT = NCBI_REFSEQ_COMP;

    const std::string INTERPRO_DEFAULT = "pfam";

    const std::string GO_DB_PATH = "/bin/go_term.entp";
    const std::string GO_TERM_FILE = "/databases/term.txt";
    const std::string TAX_SCRIPT_PATH = "/download_tax.pl";
    const std::string TAX_DATABASE_PATH = "/databases/ncbi_tax.entp";
    const std::string TAX_BIN_PATH = "/bin/ncbi_tax_bin.entp";
    const std::string BIN_PATH = "bin/";
    const std::string NCBI_INDEX_PATH = "/bin/ncbi_";

    const std::string DIAMOND_PATH_EXE = "/libs/diamond-0.8.31/bin/diamond";
    const std::string DIAMOND_INDX_OUT_PATH = "outfiles/diamond/diamond_index.out";
    const std::string DIAMOND_RUN_OUT_PATH = "outfiles/diamond/blastx_";
    const std::string SIM_SEARCH_OUT_PATH = "similarity_search/";
}

namespace ENTAP_EXECUTE {
    const std::string RSEM_EXE_PATH = "/libs/RSEM-1.3.0/";
    const std::string RSEM_OUT_DIR = "expression/";
    const float RSEM_FPKM_DEFAULT = 0.5;
    const int RSEM_COL_NUM = 7;
    const std::string OUTFILE_DEFAULT = "outfiles";
    const std::string FIGURE_DIR      = "figures/";


    //--------------------Frame Selection----------------------//

    const std::string GENEMARK_EXE_PATH = "/libs/gmst_linux_64/gmst.pl";

    //------------------------Ontology-------------------------//

    const std::string EGGNOG_EMAPPER_EXE = "/libs/eggnog-mapper/emapper.py";
    const std::string EGGNOG_INIT_EXE = "";
    const short EGGNOG_COL_NUM = 12;
    const std::string GO_BIOLOGICAL_FLAG = "biological_process";
    const std::string GO_CELLULAR_FLAG = "cellular_component";
    const std::string GO_MOLECULAR_FLAG = "molecular_function";
    const short EGGNOG_INT_FLAG = 0;

    const short INTERPRO_INT_FLAG = 1;
    const short INTERPRO_COL_NUM = 15;

    const std::string INTERPRO_EXE = "/libs/interproscan-5.22-61.0/interproscan.sh";

    const std::string ENTAP_OUTPUT = "entap_out/";

    const int diamond_col_num = 14;

    const std::list<std::string> INFORMATIVENESS {
            "conserved",
            "predicted",
            "unnamed",
            "hypothetical",
            "putative",
            "unidentified",
            "uncharacterized",
            "unknown",
            "uncultured",
            "uninformative"
    };

    //--------------------Graphing----------------------//

    const std::string GRAPH_FILEPATH   =        "/src/entap_graphing.py";

    const signed char GRAPH_INIT       =        -1;
    const signed char GRAPH_EXPRESSION =         0;
    const signed char GRAPH_FRAME      =         1;
    const signed char GRAPH_SIM        =         2;
    const signed char GRAPH_ONTOLOGY   =         3;


}

namespace ENTAP_STATS {
    const std::string SOFTWARE_BREAK = "----------------------------------------------\n";
}



#endif //ENTAP_ERRORFLAGS_H
