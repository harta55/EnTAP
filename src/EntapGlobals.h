/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
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
#include <boost/serialization/access.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <unordered_map>
#include "common.h"

class QuerySequence;
struct TaxEntry;
struct GoEntry;
namespace boostFS = boost::filesystem;
namespace boostPO = boost::program_options;
namespace boostAR = boost::archive;

//**************************************************************

//******************** Defines/Macros **************************

#define PATHS(x,y)      (boostFS::path(x) / boostFS::path(y)).string()
#define NCBI_UNIPROT    0       // Compiler flag for future feature
#define DEBUG           1
#define FASTA_FLAG      ">"

//**************************************************************


//***************** Global Prototype Functions *****************
int execute_cmd(std::string,std::string);
int execute_cmd(std::string);
std::string generate_command(std::unordered_map<std::string,std::string>&,
                             std::string);
std::string float_to_string(fp64);
std::string float_to_sci(fp64, int);
//**************************************************************


//**************** Global Structures/Typedefs ******************
typedef std::unordered_map<std::string, QuerySequence*> QUERY_MAP_T;
typedef std::unordered_map<std::string, TaxEntry> tax_serial_map_t;
typedef std::unordered_map<std::string, GoEntry> go_serial_map_t;
typedef std::map<std::string,std::vector<std::string>> go_format_t;
typedef std::vector<std::string> databases_t;   // Standard database container

typedef std::pair<std::string,int> count_pair;
struct compair {
    bool operator ()(count_pair const& one, count_pair const& two) const {
        return one.second > two.second;
    }
};

enum ExecuteStates {
    INIT = 0,
    EXPRESSION_FILTERING,
    FRAME_SELECTION,
    FILTER,
    SIMILARITY_SEARCH,
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
extern std::string TAX_DB_PATH_TEXT;// Text
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

    const uint16 FRAME_FLAG_GENEMARK = 0;
    const uint16 EXP_FLAG_RSEM       = 0;
    const uint16 SIM_SEARCH_FLAG_DIAMOND = 0;

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

namespace UInput {
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
    extern const std::string INPUT_FLAG_SINGLE_END;
    extern const std::string INPUT_FLAG_THREADS;
    extern const std::string INPUT_FLAG_UNINFORM;
    extern const std::string INPUT_FLAG_NOCHECK;
}

namespace ENTAP_STATS {
    const std::string SOFTWARE_BREAK = "------------------------------------------------------\n";
}

#endif //ENTAPGLOBALS_H
