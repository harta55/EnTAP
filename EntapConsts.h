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
    const int E_CONFIG_PARSE = 12;
    const int E_CONFIG_CREATE = 13;
    const int E_INIT_TAX_DOWN = 20;
    const int E_INIT_TAX_INDEX = 21;
    const int E_INIT_TAX_SERIAL = 22;
    const int E_INIT_INDX_DATA_NOT_FOUND = 30;
    const int E_INIT_INDX_DATABASE = 31;
    const int E_INIT_DOWNLOAD = 23;


    const int E_INIT_TAX_READ = 55;

    const int E_RUN_EXECUTION_PATHS         = 105;
    const int E_RUN_VERIFY_DATABASES        = 106;
    const int E_RUN_GENEMARK                = 100;
    const int E_RUN_RSEM_VALIDATE           = 110;
    const int E_RUN_RSEM_CONVERT            = 111;
    const int E_RUN_RSEM_EXPRESSION         = 112;
    const int E_RUN_FILTER                  = 120;
}

namespace ENTAP_CONFIG {

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


    const std::string TAX_SCRIPT_PATH = "/download_tax.pl";
    const std::string TAX_DATABASE_PATH = "/databases/ncbi_tax.entp";
    const std::string TAX_BIN_PATH = "/bin/ncbi_tax_bin.entp";
    const std::string BIN_PATH = "bin/";


    const double E_VALUE = 1e-5;

    const std::string DEBUG_FILENAME = "debug.txt";
    const std::string LOG_FILENAME = "log_file.txt";
    //------------------USER INPUTS-----------------------//
    const std::string INPUT_FLAG_ALIGN = "align";
    const std::string INPUT_FLAG_RUNPROTEIN = "runP";
    const std::string INPUT_FLAG_RUNNUCLEOTIDE = "runN";
    const std::string INPUT_FLAG_OVERWRITE = "overwrite";
    const std::string INPUT_FLAG_NCBI_1 = "ncbi";
    const std::string INPUT_FLAG_NCBI_2 = "N";
    const std::string INPUT_FLAG_UNIPROT = "uniprot";


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

    const std::string NCBI_INDEX_PATH = "/bin/ncbi_";

    //------------------FTP PATHS-----------------------//
    const std::string UNIPROT_FTP_SWISS = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";
    const std::string UNIPROT_FTP_TREMBL = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz";


    const std::string DIAMOND_PATH_EXE = "/libs/diamond-0.8.31/bin/diamond";
    const std::string DIAMOND_INDX_OUT_PATH = "outfiles/diamond/diamond_index.out";
    const std::string DIAMOND_RUN_OUT_PATH = "outfiles/diamond/blastx_";
    const std::string DIAMOND_DIR = "diamond/";
}

namespace ENTAP_EXECUTE {
    const std::string RSEM_EXE_PATH = "/libs/RSEM-1.3.0/";
    const std::string RSEM_OUT_DIR = "rsem/";
    const float RSEM_FPKM_DEFAULT = 0.5;
    const int RSEM_COL_NUM = 7;

    const std::string OUTFILE_DEFAULT = "outfiles";
    const std::string GENEMARK_EXE_PATH = "/libs/gmst_linux_64/gmst.pl";
    const std::string GENEMARK_LOG_FILE = "gms.log";
    const std::string GENEMARK_OUT_PATH = "genemark/";
    const std::string GENEMARK_OUT_UNSELECTED = "genemark_frames_not_found.fasta";

    const std::string ENTAP_OUTPUT = "entap_out/";

    const int diamond_col_num = 13;
    const int diamond_e_col = 10;
}



#endif //ENTAP_ERRORFLAGS_H
