/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
 *
 * This file is part of EnTAP.
 *
 * EnTAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EnTAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef ENTAPGLOBALS_H
#define ENTAPGLOBALS_H

//*********************** Includes *****************************
#include <vector>
#include <list>
#include <boost/serialization/access.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <unordered_map>
#include "common.h"

namespace boostFS = boost::filesystem;
namespace boostPO = boost::program_options;
namespace boostAR = boost::archive;
class QuerySequence;

//**************************************************************

//*********************** Defines ******************************

#define PATHS(x,y)      (boostFS::path(x) / boostFS::path(y)).string()
#define NCBI_UNIPROT    0       // Compiler flag for future feature
#define DEBUG           1
#define FILE_APPEND     std::ios::out | std::ios::app
#define FASTA_FLAG      ">"

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

typedef std::map<std::string, QuerySequence> QUERY_MAP_T;

typedef std::pair<std::string,int> count_pair;
struct compair {
    bool operator ()(count_pair const& one, count_pair const& two) const {
        return one.second > two.second;
    }
};

enum ExecuteStates {
    INIT,
    RSEM,
    FRAME_SELECTION,
    FILTER,
    DIAMOND_RUN,
    DIAMOND_PARSE,
    GENE_ONTOLOGY,
    EXIT
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
extern std::string TAX_DB_PATH;     // binary
extern std::string TAX_DOWNLOAD_EXE;
extern std::string GO_DB_PATH;      // binary
extern std::string GRAPHING_EXE;


namespace ENTAP_EXECUTE {
    //------------------------Ontology-------------------------//
    extern const std::string GO_BIOLOGICAL_FLAG ;
    extern const std::string GO_CELLULAR_FLAG;
    extern const std::string GO_MOLECULAR_FLAG;
    const uint16 ONTOLOGY_MIN      = 0;
    const uint16 EGGNOG_INT_FLAG   = 0;
    const uint16 INTERPRO_INT_FLAG = 1;
    const uint16 ONTOLOGY_MAX      = 1;

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
    extern const std::string HEADER_INFORM;

    // EggNOG Header Information
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

    // Interpro Header Information
    extern const std::string HEADER_INTER_GO_BIO;
    extern const std::string HEADER_INTER_GO_CELL;
    extern const std::string HEADER_INTER_GO_MOLE;
    extern const std::string HEADER_INTER_PATHWAY;
    extern const std::string HEADER_INTER_INTERPRO;
    extern const std::string HEADER_INTER_DATA_TYPE;
    extern const std::string HEADER_INTER_DATA_TERM;
    extern const std::string HEADER_INTER_EVAL;
}


namespace ENTAP_CONFIG {

    extern const std::string ENTAP_VERSION ;
    extern const std::string DEBUG_FILENAME;
    extern const std::string LOG_FILENAME  ;
    extern const std::string LOG_EXTENSION;

    //------------------USER INPUTS-----------------------//
    extern const std::string INPUT_FLAG_TAG;
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
    extern const std::string INPUT_FLAG_GRAPH;
    extern const std::string INPUT_FLAG_TRIM;
    extern const std::string INPUT_FLAG_STATE;
    extern const std::string INPUT_FLAG_PAIRED_END;

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
    extern const std::string TAX_DB_DEFAULT    ;
    extern const std::string BIN_PATH          ;
    extern const std::string DATABASE_DIR      ;
    extern const std::string NCBI_INDEX_PATH   ;
}

//**************************************************************


//******************* Global Constants *************************
namespace ENTAP_ERR {
    const uint16 E_INPUT_PARSE                 = 10;
    const uint16 E_SUCCESS                     = 0;
    const uint16 E_CONFIG_PARSE                = 12;
    const uint16 E_CONFIG_CREATE               = 13;
    const uint16 E_CONFIG_CREATE_SUCCESS       = 14;
    const uint16 E_INIT_TAX_DOWN               = 20;
    const uint16 E_INIT_TAX_INDEX              = 21;
    const uint16 E_INIT_TAX_SERIAL             = 22;
    const uint16 E_INIT_INDX_DATA_NOT_FOUND    = 30;
    const uint16 E_INIT_INDX_DATABASE          = 31;
    const uint16 E_INIT_DOWNLOAD               = 23;
    const uint16 E_INIT_EGGNOG                 = 40;

    const uint16 E_INIT_TAX_READ               = 55;
    const uint16 E_INIT_GO_DOWNLOAD            = 60;
    const uint16 E_INIT_GO_UNZIP               = 61;
    const uint16 E_INIT_GO_PARSE               = 62;
    const uint16 E_INIT_GO_INDEX               = 63;

    const uint16 E_RUN_EXECUTION_PATHS         = 105;
    const uint16 E_RUN_VERIFY_DATABASES        = 106;
    const uint16 E_RUN_GENEMARK                = 100;
    const uint16 E_RUN_GENEMARK_PARSE          = 101;
    const uint16 E_RUN_GENEMARK_STATS          = 102;
    const uint16 E_RUN_GENEMARK_MOVE           = 103;
    const uint16 E_RUN_RSEM_VALIDATE           = 110;
    const uint16 E_RUN_RSEM_CONVERT            = 111;
    const uint16 E_RUN_RSEM_EXPRESSION         = 112;
    const uint16 E_RUN_RSEM_EXPRESSION_PARSE   = 113;
    const uint16 E_RUN_FILTER                  = 120;
    const uint16 E_RUN_SIM_SEARCH_FILTER       = 140;
    const uint16 E_RUN_SIM_SEARCH_RUN          = 141;
    const uint16 E_RUN_ANNOTATION              = 150;
    const uint16 E_GO_READ                     = 151;
    const uint16 E_RUN_EGGNOG                  = 160;
    const uint16 E_DATABASE_QUERY              = 161;
    const uint16 E_PARSE_EGGNOG                = 162;
    const uint16 E_RUN_INTERPRO                = 170;
    const uint16 E_PARSE_INTERPRO              = 171;
}

namespace ENTAP_STATS {
    const std::string SOFTWARE_BREAK = "----------------------------------------------\n";
}

//***********************************************




#endif //ENTAPGLOBALS_H
