//
// Created by harta55 on 2/1/17.
//

#ifndef ENTAP_ERRORFLAGS_H
#define ENTAP_ERRORFLAGS_H
#include <vector>

// TODO temporary, change
namespace ENTAP_ERR {
    const int E_INPUT_PARSE = 10;
    const int E_SUCCESS = 11;
    const int E_INIT_TAX_DOWN = 20;
    const int E_INIT_TAX_INDEX = 21;
    const int E_INIT_TAX_SERIAL = 22;
    const int E_INIT_INDX_DATA_NOT_FOUND = 30;


    const int E_INIT_TAX_READ = 55;

    const int E_RUN_GENEMARK                = 100;
    const int E_RUN_RSEM_VALIDATE           = 110;
    const int E_RUN_RSEM_CONVERT            = 111;
    const int E_RUN_RSEM_EXPRESSION         = 112;
}

namespace ENTAP_CONFIG {
    const std::string TAX_SCRIPT_PATH = "/download_tax.pl";
    const std::string TAX_DATABASE_PATH = "/databases/ncbi_tax.entp";
    const std::string TAX_BIN_PATH = "/bin/ncbi_tax_bin.entp";
    const std::string BIN_PATH = "bin/";


    const double E_VALUE = 1e-5;
    //------------------USER INPUTS-----------------------//
    const std::string INPUT_UNIPROT_SWISS = "swiss";
    const std::string INPUT_UNIPROT_UR100 = "ur100";
    const std::string INPUT_UNIPROT_UR90 = "ur90";
    const std::string INPUT_UNIPROT_TREMBL = "trembl";
    const std::string INPUT_UNIPROT_NULL = "null";
    const std::string INPUT_UNIPROT_DEFAULT = INPUT_UNIPROT_SWISS;

    const std::string UNIPROT_FTP_SWISS = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";
    const std::string UNIPROT_BASE_PATH = "/databases/uniprot_";
    const std::string UNIPROT_INDEX_PATH = "/bin/uniprot_";

    const std::string NCBI_NONREDUNDANT = "nr";
    const std::string NCBI_BASE_PATH = "/databases/ncbi_";
    const std::string NCBI_REFSEQ_COMP = "refseq-c";
    const std::string NCBI_REFSEQ_PLANT = "refseq-p";
    const std::string NCBI_NULL = "null";
    const std::string NCBI_DEFAULT = NCBI_REFSEQ_COMP;

    const std::string NCBI_INDEX_PATH = "/bin/ncbi_";

    const std::string DIAMOND_PATH_EXE = "../libs/diamond-0.8.31/bin/diamond";
    const std::string DIAMOND_INDX_OUT_PATH = "outfiles/diamond/diamond_index.out";
    const std::string DIAMOND_RUN_OUT_PATH = "outfiles/diamond/blastx_";
    const std::string DIAMOND_DIR = "outfiles/diamond";

}

namespace ENTAP_EXECUTE {
    const std::string RSEM_EXE_PATH = "/libs/RSEM-1.3.0/";
    const float RSEM_FPKM_DEFAULT = 0.5;
    const int RSEM_COL_NUM = 7;


    const std::string EXECUTE_OUT_PATH = "outfiles/";
    const std::string GENEMARK_EXE_PATH = "/libs/gmst_linux_64/gmst.pl";

    const int diamond_col_num = 13;
    const int diamond_e_col = 10;
}



#endif //ENTAP_ERRORFLAGS_H
