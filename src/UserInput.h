/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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

#ifndef ENTAP_USERINPUT_H
#define ENTAP_USERINPUT_H

//*********************** Includes *****************************
#include "EntapGlobals.h"
#include "FileSystem.h"
#include "config.h"
#include "database/EntapDatabase.h"

#ifdef USE_BOOST
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#else
#include <tclap/CmdLine.h>
#include <boost/any.hpp>          // Include any boost library with tclap
#endif
//**************************************************************

// Define User Input types whenever an input is called MUST be one of these
// boost::any refers to these types
typedef vect_str_t ent_input_multi_str_t;
typedef vect_fp64_t ent_input_multi_fp_t;
typedef vect_uint16_t ent_input_multi_int_t;
typedef std::string   ent_input_str_t;
typedef fp64 ent_input_fp_t;
typedef uint16 ent_input_uint_t;

// WARNING order must match ENTAP_INPUT_FLAGS[] in UserInput.cpp
typedef enum {
    INPUT_FLAG_UNUSED=0,

    /* General Commands */
    INPUT_FLAG_OUTPUT_DIR,
    INPUT_FLAG_CONFIG,
    INPUT_FLAG_RUNPROTEIN,
    INPUT_FLAG_RUNNUCLEOTIDE,
    INPUT_FLAG_OVERWRITE,
    INPUT_FLAG_INI_FILE,
//    INPUT_FLAG_HELP,      // Native to TCLAP
//    INPUT_FLAG_VERSION,   // Native to TCLAP
    INPUT_FLAG_TRANSCRIPTOME,
    INPUT_FLAG_DATABASE,
    INPUT_FLAG_GRAPH,
    INPUT_FLAG_NO_TRIM,
    INPUT_FLAG_THREADS,
    INPUT_FLAG_STATE,
    INPUT_FLAG_NOCHECK,
    INPUT_FLAG_OUTPUT_FORMAT,

    /* EnTAP Commands */
    INPUT_FLAG_ENTAP_DB_BIN,
    INPUT_FLAG_ENTAP_DB_SQL,
    INPUT_FLAG_ENTAP_GRAPH,
    INPUT_FLAG_ENTAP_HEADERS,

    /* Configuration Commands */
    INPUT_FLAG_DATABASE_GENERATE,
    INPUT_FLAG_DATABASE_TYPE,

    /* Expression Analysis Commands */
    INPUT_FLAG_FPKM,
    INPUT_FLAG_ALIGN,
    INPUT_FLAG_SINGLE_END,

    /* Expression ANalysis - RSEM Commands */
    INPUT_FLAG_RSEM_CALC_EXPRES,
    INPUT_FLAG_RSEM_SAM_VALID,
    INPUT_FLAG_RSEM_PREP_REF,
    INPUT_FLAG_RSEM_CONVERT_SAM,

    /* Frame Selection Commands */
    INPUT_FLAG_COMPLETE,
    INPUT_FLAG_FRAME_SELECTION,

    /* Frame Selection - GeneMarkST Commands */
    INPUT_FLAG_GENEMARKST_EXE,

    /* Frame Selection - TransDecoder Commands */
    INPUT_FLAG_TRANS_LONGORF_EXE,
    INPUT_FLAG_TRANS_PREDICT_EXE,
    INPUT_FLAG_TRANS_MIN_PROTEIN,
    INPUT_FLAG_TRANS_NO_REFINE_STARTS,

    /* Similarity Search Commands */
    INPUT_FLAG_DIAMOND_EXE,
    INPUT_FLAG_SPECIES,
    INPUT_FLAG_QCOVERAGE,
    INPUT_FLAG_TCOVERAGE,
    INPUT_FLAG_CONTAMINANT,
    INPUT_FLAG_E_VALUE,
    INPUT_FLAG_UNINFORMATIVE,

    /* Ontology Commands */
    INPUT_FLAG_ONTOLOGY,
    INPUT_FLAG_GO_LEVELS,

    /* Ontology Commands - EggNOG */
    INPUT_FLAG_EGG_SQL_DB,
    INPUT_FLAG_EGG_DMND_DB,

    /* Ontology Commands - InterProScan */
    INPUT_FLAG_INTERPRO_EXE,
    INPUT_FLAG_INTERPRO,

    /* Ontology Commands - BUSCO */
    INPUT_FLAG_BUSCO_EXE,
    INPUT_FLAG_BUSCO_DATABASE,
    INPUT_FLAG_BUSCO_EVAL,


    INPUT_FLAG_MAX

} ENTAP_INPUT_FLAGS;


class UserInput {

public:
    UserInput(int argc, const char** argv, FileSystem*fileSystem);
    ~UserInput();


    static constexpr uint16 MAX_GO_LEVEL = 12;
    static constexpr uint16 MIN_GO_LEVEL = 0;

    static std::string getBIN_PATH_DEFAULT();
    static const std::string &getENTAP_DATABASE_BIN_DEFAULT();
    static const std::string &getENTAP_DATABASE_SQL_DEFAULT();
    static const std::string &getEGG_SQL_DB_FILENAME();
    static const std::string &getEGG_DMND_FILENAME();
    static const std::string &getEGG_SQL_DB_DEFAULT();
    static const std::string &getEGG_DMND_DEFAULT();
    static std::string getDATABASE_DIR_DEFAULT();

    bool has_input(ENTAP_INPUT_FLAGS input);
    void parse_ini(std::string &ini_path);
    bool verify_user_input();
    int get_supported_threads();
    std::queue<char> get_state_queue();
    std::string get_target_species_str();
    vect_str_t get_contaminants();
    vect_str_t get_uninformative_vect();
    std::string get_user_transc_basename();
    ent_input_str_t get_entap_database_path(EntapDatabase::DATABASE_TYPE type);
    std::vector<FileSystem::ENT_FILE_TYPES> get_user_output_types();
    bool run_frame_selection(QueryData *queryData, bool &run_frame_selection);
    bool run_expression_filtering();


    template<class T>
    T get_user_input(ENTAP_INPUT_FLAGS key) {
     // Use TCLAP
        if (has_input(key)) {
            return boost::any_cast<T>(mUserInputs[key].parsed_value);
        } else {
            return T();
        }
    }

private:

    typedef enum {
        ENT_COMMAND_LINE,
        ENT_INI_FILE,
        ENT_INPUT_FUTURE    // Flag for a future feature, will not be in other input methods yet
    } ENT_INPUT_TYPES;

    // WARNING order must match VAR_TYPE_STR
    typedef enum {
        ENT_INI_VAR_STRING=0,
        ENT_INI_VAR_MULTI_STRING,
        ENT_INI_VAR_INT,
        ENT_INI_VAR_MULTI_INT,
        ENT_INI_VAR_FLOAT,
        ENT_INI_VAR_MULTI_FLOAT,
        ENT_INI_VAR_BOOL,

        ENT_INI_MAX_TYPES

    } ENT_INI_VAR_TYPES;

    // WARNING order must match ENT_INI_VAR_TYPES
    const std::string VAR_TYPE_STR[ENT_INI_MAX_TYPES] {
        "string",
        "list (string)",
        "integer",
        "list (integer)",
        "decimal",
        "list (decimal)",
        "boolean (true/false)"
    };


    struct EntapINIEntry {
        std::string category;
        std::string input;
        std::string short_input;
        std::string description;
        std::string example;
        ENT_INI_VAR_TYPES   var_type;
        boost::any          default_value;
        ENT_INPUT_TYPES     input_type;
        boost::any          parsed_value;
    };

    typedef enum {
        SPECIES,
        CONTAMINANT
    } SPECIES_FLAGS;

    void parse_future_inputs();
    void parse_arguments_tclap(int, const char **);
    void print_user_input();
    EntapINIEntry* check_ini_key(std::string &);
    void verify_databases(bool);
    void verify_species (SPECIES_FLAGS, EntapDatabase*);
    void process_user_species(std::string&);
    void verify_software_paths(std::string &state, bool is_protein, bool is_execution, QueryData *pQuery_data);
    void generate_ini_file(std::string& ini_path);

    const uint16 MAX_GO_LEVELS_SELECTED        = 5; // Max number of GO levels user can select to output
    const uint16 TRANS_MIN_PROTEIN_MIN         = 0;
    const fp32 COVERAGE_MIN                    = 0.0;
    const fp32 COVERAGE_MAX                    = 100.0;
    const fp32 FPKM_MIN                        = 0.0;
    const fp32 FPKM_MAX                        = 100.0;
    const uint8 MAX_DATABASE_SIZE              = 5;         // Maximum number of databases allowed from user

    static const std::string RSEM_DEFAULT_EXE_DIR                  ;
    static const std::string RSEM_SAM_VALID    ;
    static const std::string RSEM_PREP_REF_EXE ;
    static const std::string RSEM_CALC_EXP_EXE ;
    static const std::string RSEM_CONV_SAM     ;
    static const std::string DEFAULT_RSEM_SAM_VALID;
    static const std::string DEFAULT_RSEM_PREP_REF;
    static const std::string DEFAULT_RSEM_CALC_EXP;
    static const std::string DEFAULT_RSEM_CONV_SAM;
    static const std::string GENEMARK_DEFAULT_EXE              ;
    static const std::string TRANSDECODER_LONG_DEFAULT_EXE     ;
    static const std::string TRANSDECODER_PREDICT_DEFAULT_EXE ;
    static const std::string DIAMOND_DEFAULT_EXE               ;
    static const std::string EGG_SQL_DB_FILENAME               ;
    static const std::string EGG_DMND_FILENAME                 ;
    static const std::string INTERPRO_DEF_EXE                  ;
    static const std::string GRAPH_SCRIPT_DEF                  ;
    static const std::string DEFAULT_BUSCO_EXE;
    static const std::string BIN_PATH_DEFAULT                  ;
    static const std::string DATABASE_DIR_DEFAULT              ;
    static const std::string ENTAP_DATABASE_SQL_FILENAME       ;
    static const std::string ENTAP_DATABASE_SERIAL_FILENAME    ;
    static const std::string ENTAP_DATABASE_BIN_DEFAULT        ;
    static const std::string ENTAP_DATABASE_SQL_DEFAULT        ;
    static const std::string EGG_SQL_DB_DEFAULT                ;
    static const std::string EGG_DMND_DEFAULT                  ;
    static const std::string DEFAULT_ENTAP_DB_BIN_INI;
    static const std::string DEFAULT_ENTAP_DB_SQL_INI;
    static const std::string DEFAULT_EGG_SQL_DB_INI;
    static const std::string DEFAULT_EGG_DMND_DB_INI;
    static const std::string DEFAULT_ENTAP_GRAPH_INI;

    static const fp64          DEFAULT_E_VALUE;
    static const fp64          DEFAULT_BUSCO_E_VALUE;
    static const uint16        DEFAULT_THREADS;
    static const fp64          RSEM_FPKM_DEFAULT;
    static const fp64          DEFAULT_QCOVERAGE;
    static const fp64          DEFAULT_TCOVERAGE;
    static const vect_str_t    DEFAULT_UNINFORMATIVE;
    static const vect_uint16_t DEFAULT_DATA_TYPE;
    static const std::string   DEFAULT_OUT_DIR;
    static const std::string   DEFAULT_STATE;
    static const uint16        DEFAULT_FRAME_SELECTION;
    static const vect_uint16_t DEFAULT_ONT_LEVELS;
    static const vect_uint16_t DEFAULT_ONTOLOGY;
    static const vect_uint16_t DEFAULT_OUT_FORMAT;
    static const uint16        DEFAULT_TRANSDECODER_MIN_PROTEIN;
    static const std::string   DEFAULT_INI_PATH;

    static const std::string ENTAP_INI_FILENAME;
    const std::string INI_FILE_BOOL_TRUE = "true";
    const std::string INI_FILE_BOOL_FALSE= "false";
    const char INI_FILE_MULTI_DELIM = ',';
    const char INI_FILE_COMMENT   = '#';
    const char INI_FILE_ASSIGN    = '=';        // Assignment char for INI file

    FileSystem *mpFileSystem;
    std::string mIniFilePath;
    bool        mIsConfig;
    static EntapINIEntry mUserInputs[];
};





#endif //ENTAP_USERINPUT_H
