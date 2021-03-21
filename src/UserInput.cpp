/******************************************************************
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
 *******************************************************************/


//*********************** Includes *****************************
#include "UserInput.h"
#include "ExceptionHandler.h"
#include "GraphingManager.h"
#include "common.h"
#include "ontology/ModEggnog.h"
#include "ontology/ModInterpro.h"
#include "version.h"
#include "database/EntapDatabase.h"
#include "ontology/ModEggnogDMND.h"
#include "config.h"
#include "similarity_search/ModDiamond.h"
#include "FrameSelection.h"
#include "frame_selection/ModTransdecoder.h"
#include "database/BuscoDatabase.h"
#include "ontology/ModBUSCO.h"
#include "frame_selection/ModGeneMarkST.h"
//**************************************************************

//*********************** Defines ******************************


/* ----------------- General Input Commands ------------------*/
#define CMD_HELP            "help"
#define CMD_SHORT_HELP      "h"
#define DESC_HELP           "Print all the help options for this version of EnTAP!"
#define CMD_CONFIG          "config"
#define DESC_CONFIG         "Configure EnTAP for execution later. If this is your first time running EnTAP run this first!" \
                            "\nThis will perform the following:\n"                           \
                            "    - Downloading EnTAP/NCBI taxonomic database\n"              \
                            "    - Downloading Gene Ontology term database\n"                \
                            "    - Formatting any database you would like for diamond\n"                                      \
                            "    - Downloading UniProt Swiss-Prot information"
#define CMD_RUN_PROTEIN     "runP"
#define DESC_RUN_PROTEIN    "Execute EnTAP functionality with 'blastp'. This means that EnTAP will use protein sequences"   \
                            " for all annotation stages. If you input a nucleotide transcriptome, they will be frame"       \
                            " selected and the subsequent protein file will be used for annotation.\n"                      \
                            "This is typically how EnTAP would be ran."
#define CMD_RUN_NUCLEO      "runN"
#define DESC_RUN_NUCLEO     "Execute EnTAP functionality with 'blastx'. This means that EnTAP will use nucleotide sequences" \
                            " for all annotation stages. If you input a nucleotide trancsriptome, that will be used to "     \
                            " annotate. Due to this, you will not be able to select this option and input a protein transcriptome."

#define DESC_INPUT_TRAN     "Path to the input transcriptome file"
#define CMD_INPUT_TRAN      "input"
#define CMD_SHORT_INPUT_TRAN "i"

#define DESC_OUTPUT_DIR     "Specify the output directory you would like the data produced by EnTAP to be saved to."
#define CMD_OUTPUT_DIR      "out-dir"

#define CMD_OVERWRITE       "overwrite"
#define DESC_OVERWRITE      "Select this option if you would like to overwrite files from a previous execution of EnTAP."   \
                            " This will DISABLE 'picking up where you left off' which enables you to continue an annotation" \
                            " from where you left off before. Refer to the documentation for more information."
#define CMD_INI_FILE        "ini"
#define DESC_INI_FILE      "[REQUIRED] Specify path to the entap_config.ini file that will be used to find all of the configuration data."

#define DESC_VERSION       "Print the version of EnTAP software you are running."
#define CMD_VERSION        "version"
#define CMD_SHORT_VERSION  "v"
#define DESC_GRAPHING       "Check whether or not your system supports graphing. This option does not require any other flags and will"     \
                            " exit EnTAP after it determined that the proper Python libraries are present"
#define CMD_GRAPHING        "graph"

#define DESC_NO_TRIM        "By default, EnTAP will trim the sequence ID to the nearest space to help with compatibility" \
                            " across software. This command will instead remove the spaces in a sequence ID rather than trimming."
#define EX_NO_TRIM          "'>TRINITY_231.1 Protein Information' will become...\n'>TRINITY_231.1ProteinInformation' \n"
#define CMD_NO_TRIM         "no-trim"

#define DESC_THREADS        "Specify the number of threads that will be used throughout EnTAP execution\n"
#define CMD_THREADS         "threads"
#define CMD_SHORT_THREADS   "t"
const uint16 UserInput::DEFAULT_THREADS                 = 1;

#define CMD_STATE           "state"
#define DESC_STATE          "Specify the state of execution (EXPERIMENTAL). More information is available in the documentation." \
                            " This flag may have undesired affects and may not run properly!"

#define DESC_NOCHECK        "Use this flag if you don't want your input to EnTAP verifed. This is not advised to use! Your run " \
                            "may fail later on if inputs are not checked."
#define CMD_NOCHECK         "no-check"

#define CMD_OUTPUT_FORMAT   "output-format"
#define DESC_OUTPUT_FORMAT  "Specify the output format for the processed alignments."   \
                            "Multiple flags can be specified:\n"                          \
                            "    1. TSV Format (default)\n"                             \
                            "    2. CSV Format\n"                                       \
                            "    3. FASTA Amino Acid (default)\n"                       \
                            "    4. FASTA Nucleotide (default)\n"                       \
                            "    5. Gene Enrichment Sequence ID vs. Effective Length TSV (default)\n"\
                            "    6. Gene Enrichment Sequence ID vs. GO Term TSV (default)"
const vect_uint16_t UserInput::DEFAULT_OUT_FORMAT       =vect_uint16_t{FileSystem::ENT_FILE_DELIM_TSV,
                                                                       FileSystem::ENT_FILE_FASTA_FAA,
                                                                       FileSystem::ENT_FILE_FASTA_FNN,
                                                                       FileSystem::ENT_FILE_GENE_ENRICH_EFF_LEN,
                                                                       FileSystem::ENT_FILE_GENE_ENRICH_GO_TERM};

/* -------------------- EnTAP Commands -----------------------*/
#define CMD_ENTAP_DB_BIN    "entap-db-bin"
#define DESC_ENTAP_DB_BIN   "Path to the EnTAP binary database"
#define CMD_ENTAP_DB_SQL    "entap-db-sql"
#define DESC_ENTAP_DB_SQL   "Path to the EnTAP SQL database (not needed if you are using the binary database)"
#define CMD_ENTAP_GRAPH_PATH "entap-graph"
#define DESC_ENTAP_GRAPH_PATH "Path to the EnTAP graphing script (entap_graphing.py)"

/* ---------------- Configuration Commands -------------------*/
#define CMD_DATA_GENERATE   "data-generate"
#define DESC_DATA_GENERATE  "Specify whether you would like to generate EnTAP databases locally instead of downloading them." \
                            " By default, EnTAP will download the databases. This may be used if you are experiencing "\
                            " errors with the default process."

#define CMD_DATABASE_TYPE   "data-type"
#define DESC_DATABASE_TYPE  "Specify which EnTAP database you would like to download/generate or use throughout execution." \
                            " Only one is required.\n"     \
                            "    0. Serialized Database (default)\n"                    \
                            "    1. SQLITE Database\n"                                  \
                            "It is advised to use the default Serialized Database as this is fastest."

/* ---------------- Expression Analysis Commands -------------*/
#define CMD_FPKM            "fpkm"
#define DESC_FPKM           "Specify the FPKM threshold with expression analysis. EnTAP will filter out transcripts below" \
                            " this value. (default: 0.5)"
const fp64   UserInput::RSEM_FPKM_DEFAULT               = 0.5;

#define DESC_ALIGN_FILE     "Specify the path to the BAM/SAM file for expression analysis"
#define CMD_ALIGN_FILE      "align"
#define CMD_SHORT_ALIGN_FILE "a"
#define CMD_SINGLE_END      "single-end"
#define DESC_SINGLE_END     "Specify this flag if your BAM/SAM file was generated "    \
                            "through single-end reads\n"                                \
                            "Note: this is only required in expression analysis\n"      \
                            "Default: paired-end"
#define DESC_RSEM_CALC_EXP  "Execution method of RSEM Calculate Expression."
#define EX_RSEM_CALC_EXP    "Example: rsem-calculate-expression"
#define CMD_RSEM_CALC_EXP   "rsem-calculate-expression"
#define DESC_RSEM_SAM_VALID "Execution method of RSEM SAM Validate."
#define EX_RSEM_SAM_VALID   "Example: rsem-sam-validator"
#define CMD_RSEM_SAM_VALID  "rsem-sam-validator"
#define DESC_RSEM_PREP_REF  "Execution method of RSEM Prep Reference."
#define EX_RSEM_PREP_REF    "Example: rsem-prepare-reference"
#define CMD_RSEM_PREP_REF   "rsem-prepare-reference"
#define DESC_RSEM_CONV_SAM  "Execution method of RSEM Convert SAM"
#define EX_RSEM_CONV_SAM    "Example: convert-sam-for-rsem"
#define CMD_RSEM_CONV_SAM   "convert-sam-for-rsem"

/* ----------------- Frame Selection Commands ----------------*/
#define CMD_GENEMARKST_EXE  "genemarkst-exe"
#define DESC_GENEMARKST_EXE "Method to execute GeneMarkST. This may be the path to the executable."
#define CMD_TRANS_LONG_EXE  "transdecoder-long-exe"
#define DESC_TRANS_LONG_EXE "Method to execute TransDecoder.LongOrfs. This may be the path to "\
                            "the executable or simply TransDecoder.LongOrfs"
#define CMD_TRANS_PREDICT_EXE "transdecoder-predict-exe"
#define DESC_TRANS_PREDICT_EXE "Method to execute TransDecoder.Predict. This may be the path to "\
                            "the executable or simply TransDecoder.Predict"
#define CMD_FRAME_SELECTION_FLAG "frame-selection"
#define DESC_FRAME_SELECTION_FLAG "Specify the Frame Selection software you would like "\
                            "to use. Only one flag can be specified.\n"                 \
                            "Specify flags as follows:\n"                               \
                            "    1. GeneMarkS-T\n"                            \
                            "    2. Transdecoder (default)"
#define DESC_COMPLETE_PROT  "Select this option if all of your sequences are complete " \
                            "proteins.\n"                                               \
                            "At this point, this option will merely flag the sequences in your output file"
#define CMD_COMPLETE_PROT   "complete"
#define DESC_TRANS_MIN_FLAG "Transdecoder only. Specify the minimum protein length"
#define CMD_TRANS_MIN_FLAG "transdecoder-m"
#define CMD_TRANS_NO_REF_START "transdecoder-no-refine-starts"
#define DESC_TRANS_NO_REF_START "Specify this flag if you would like to pipe the TransDecoder " \
                                "command '--no_refine_starts' when it is executed. Default: False\n"           \
                                "This will 'start refinement identifies potential start codons for " \
                                "5' partial ORFs using a PWM, process on by default.' "
/* ------------------ Similarity Search Commands -------------*/
#define CMD_DATABASE        "database"
#define CMD_SHORT_DATABASE  "d"
#define DESC_DATABASE       "Provide the paths to the databases you would like to use for either 'run' or 'configuration'." \
                            "\nFor running/execution:\n"                                                                    \
                            "    - Ensure the databases selected are in a DIAMOND configured format with an extension of .dmnd" \
                            "\nFor configuration:\n"                                                                             \
                            "    - Ensure the databases are in a typical FASTA format"\
                            "\nNote: if your databases do not have the typical NCBI or UniProt header format, taxonomic "   \
                            " information and filtering may not be utilized. Refer to the documentation to see how to properly format" \
                            " any data."

#define DESC_DIAMOND_EXE     "Method to execute DIAMOND. This can be a path to the executable or simply 'diamond' if installed globally."
#define CMD_DIAMOND_EXE     "diamond-exe"
#define DESC_TAXON          "Specify the type of species/taxon you are analyzing and would like alignments closer in taxonomic relevance" \
                            " to be favored (based on NCBI Taxonomic Database)\n"              \
                            "Note: replace all spaces with underscores '_'"
#define CMD_TAXON           "taxon"
#define CMD_QCOVERAGE       "qcoverage"
#define DESC_QCOVERAGE      "Select the minimum query coverage to be allowed during similarity searching"
const fp64   UserInput::DEFAULT_QCOVERAGE               = 50.0;

#define CMD_TCOVERAGE       "tcoverage"
#define DESC_TCOVERAGE      "Select the minimum target coverage to be allowed during similarity searching"
const fp64   UserInput::DEFAULT_TCOVERAGE               = 50.0;

#define CMD_CONTAMINANT     "contam"
#define CMD_SHORT_CONTAMINANT "c"
#define DESC_CONTAMINANT    "Specify the contaminants you would like to flag for similarity searching. Contaminants can" \
                            " be selected by species or through a specific taxon (insecta) from the NCBI Taxonomy Database." \
                            " If your taxon is more than one word just replace the spaces with underscores (_)." \
                            "\nNote: since hits are based upon a multitide of factors, a contaminant might end up being the best " \
                            "hit for an alignment. In this scenario, EnTAP will flag the contaminant and it can be removed"\
                            " if you would like."

#define CMD_EVAL            "e-value"
#define CMD_SHORT_EVAL      "e"
#define DESC_EVAL           "Specify the E-Value that will be used as a cutoff during similarity searching."
const fp64   UserInput::DEFAULT_E_VALUE                 = 1e-5;

#define CMD_UNINFORMATIVE   "uninformative"
#define DESC_UNINFORMATIVE  "List of keywords that should be used to specify uninformativeness of hits during similarity searching. "  \
                            "Generally something along the lines of 'hypothetical' or 'unknown' are used. Each term should be separated by a comma (,) " \
                            "This can be used if you would like to tag certain descriptions or "\
                            "would like to weigh certain alignments differently (see full documentation)"\
                            "\nExample (defaults):\n"          \
                            "conserved, predicted, unknown, hypothetical, putative, unidentified, uncultured, uninformative, unnamed"

// Enter as lowercase
const vect_str_t UserInput::DEFAULT_UNINFORMATIVE       = vect_str_t {
        "conserved",
        "predicted",
        "unknown",
        "unnamed",
        "hypothetical",
        "putative",
        "unidentified",
        "uncharacterized",
        "uncultured",
        "uninformative"
};

/* -------------------- Ontology Commands --------------------*/
#define DESC_EGGNOG_DMND     "Path to EggNOG DIAMOND configured database that was generated during the Configuration stage."
#define CMD_EGGNOG_DMND     "eggnog-dmnd"

#define DESC_EGGNOG_SQL     "Path to the EggNOG SQL database that was downloaded during the Configuration stage."
#define CMD_EGGNOG_SQL      "eggnog-sql"

#define DESC_INTERPRO_EXE   "Execution method of InterProScan. This is how InterProScan is generally ran on your system. " \
                            " It could be as simple as 'interproscan.sh' depending on if it is globally installed."
#define CMD_INTERPRO_EXE    "interproscan-exe"

#define CMD_INTER_DATA      "protein"
#define DESC_INTER_DATA     "Select which databases you would like for InterProScan. "    \
                            "Databases must be one of the following:\n"                 \
                            "    -tigrfam\n"                                            \
                            "    -sfld\n"                                               \
                            "    -prodom\n"                                             \
                            "    -hamap\n"                                              \
                            "    -pfam\n"                                               \
                            "    -smart\n"                                              \
                            "    -cdd\n"                                                \
                            "    -prositeprofiles\n"                                    \
                            "    -prositepatterns\n"                                    \
                            "    -superfamily\n"                                        \
                            "    -prints\n"                                             \
                            "    -panther\n"                                            \
                            "    -gene3d\n"                                             \
                            "    -pirsf\n"                                              \
                            "    -coils\n"                                              \
                            "    -morbidblite\n"                                        \
                            "Make sure the database is downloaded, EnTAP will not check!"
#define EX_INTER_DATA       "--" CMD_INTER_DATA " tigrfam " "--" CMD_INTER_DATA " pfam"

#define CMD_ONTOLOGY_FLAG   "ontology"
#define DESC_ONTOLOGY_FLAG  " Specify the ontology software you would like to use\n"     \
                            "Note: it is possible to specify more than one! Just use"   \
                            "multiple --ontology flags\n"                               \
                            "Specify flags as follows:\n"                               \
                            "    0. EggNOG (default)\n"                                 \
                            "    1. InterProScan"
#define CMD_GO_LEVELS      "level"
#define DESC_ONT_LEVELS     "Specify the Gene Ontology levels you would like printed\n" \
                            "A level of 0 means that every term will be printed, while a level of 1 or higher\n"       \
                            "means that that level and anything higher than it will be printed\n"   \
                            "It is possible to specify multiple flags as well\n"                     \
                            "Example/Defaults: --level 0 --level 1"

/* BUSCO */
#define CMD_BUSCO_EXE      "busco-exe"
#define DESC_BUSCO_EXE     "Specify the execution method of BUSCO."
#define EX_BUSCO_EXE       "Example: run_BUSCO.py"
#define CMD_BUSCO_DATABASE "busco-database"
#define DESC_BUSCO_DATABASE "Specify the BUSCO/OrthoDB databases you would like to download. They \n"  \
                            "can be specified by their type (eukaryota) or through a direct link to \n"\
                            "the database (http://www...)."
#define CMD_BUSCO_EVAL     "busco-eval"
#define DESC_BUSCO_EVAL    "Minimum E-Value for BUSCO related BLAST searches."
const fp64   UserInput::DEFAULT_BUSCO_E_VALUE           = 1e-5;

#define INI_FRAME_GENEMARK "frame_selection-genemarks-t"
#define INI_FRAME_TRANSDECODER "frame_selection-transdecoder"
#define INI_GENERAL "general"
#define INI_CONFIG "configuration"
#define INI_EXPRESSION "expression_analysis"
#define INI_EXP_RSEM "expression_analysis-rsem"
#define INI_FRAME   "frame_selection"
#define INI_SIM_SEARCH "similarity_search"
#define INI_ONT_INTERPRO "ontology-interproscan"
#define INI_ONTOLOGY "ontology"
#define INI_ONT_EGGNOG "ontology-eggnog"
#define INI_ENTAP "entap"
#define INI_TRANSC_BUSCO "ontology-busco"

#define ENTAP_INI_NULL_STR_VECT vect_str_t()
#define ENTAP_INI_NULL_INT_VECT vect_uint16_t()
#define ENTAP_INI_NULL_FLT_VECT vect_fp64_t()
#define ENTAP_INI_NULL_INT (0u)
#define ENTAP_INI_NULL_FLT (0.0f)
#define ENTAP_INI_NULL std::string()
#define ENTAP_INI_NULL_VAL boost::any()

//**************************************************************
const std::string UserInput::ENTAP_INI_FILENAME         = "entap_config.ini";

const vect_uint16_t UserInput::DEFAULT_DATA_TYPE        = vect_uint16_t{EntapDatabase::ENTAP_SERIALIZED};
const std::string UserInput::DEFAULT_STATE              ="+";
const uint16 UserInput::DEFAULT_FRAME_SELECTION         = FRAME_TRANSDECODER;
const vect_uint16_t UserInput::DEFAULT_ONT_LEVELS       =vect_uint16_t{0,1};
const vect_uint16_t UserInput::DEFAULT_ONTOLOGY         =vect_uint16_t{ONT_EGGNOG_DMND};

const uint16      UserInput::DEFAULT_TRANSDECODER_MIN_PROTEIN = 100;

const std::string UserInput::TRANSDECODER_LONG_DEFAULT_EXE    = "TransDecoder.LongOrfs";
const std::string UserInput::TRANSDECODER_PREDICT_DEFAULT_EXE = "TransDecoder.Predict";
const std::string UserInput::DIAMOND_DEFAULT_EXE              = PATHS(FileSystem::get_exe_dir(),"/libs/diamond-0.9.9/bin/diamond");
const std::string UserInput::EGG_SQL_DB_FILENAME              = "eggnog.db";
const std::string UserInput::EGG_DMND_FILENAME                = "eggnog_proteins.dmnd";
const std::string UserInput::INTERPRO_DEF_EXE                 = "interproscan.sh";
const std::string UserInput::GRAPH_SCRIPT_DEF                 = "/src/entap_graphing.py";
const std::string UserInput::BIN_PATH_DEFAULT                 = "/bin";
const std::string UserInput::DATABASE_DIR_DEFAULT             = "/databases";
const std::string UserInput::ENTAP_DATABASE_SQL_FILENAME      = "entap_database.db";
const std::string UserInput::ENTAP_DATABASE_SERIAL_FILENAME   = "entap_database.bin";
const std::string UserInput::ENTAP_DATABASE_BIN_DEFAULT       = PATHS(BIN_PATH_DEFAULT, ENTAP_DATABASE_SERIAL_FILENAME);
const std::string UserInput::ENTAP_DATABASE_SQL_DEFAULT       = PATHS(DATABASE_DIR_DEFAULT, ENTAP_DATABASE_SQL_FILENAME);
const std::string UserInput::EGG_SQL_DB_DEFAULT               = PATHS(DATABASE_DIR_DEFAULT, EGG_SQL_DB_FILENAME);
const std::string UserInput::EGG_DMND_DEFAULT                 = PATHS(BIN_PATH_DEFAULT, EGG_DMND_FILENAME);

// INI file path defaults using static filesystem
const std::string UserInput::RSEM_SAM_VALID        = "rsem-sam-validator";
const std::string UserInput::RSEM_PREP_REF_EXE     = "rsem-prepare-reference";
const std::string UserInput::RSEM_CALC_EXP_EXE     = "rsem-calculate-expression";
const std::string UserInput::RSEM_CONV_SAM         = "convert-sam-for-rsem";
const std::string UserInput::RSEM_DEFAULT_EXE_DIR  = PATHS(FileSystem::get_exe_dir(),"/libs/RSEM-1.3.3/");   // Directory
const std::string UserInput::DEFAULT_RSEM_SAM_VALID= PATHS(RSEM_DEFAULT_EXE_DIR, RSEM_SAM_VALID);
const std::string UserInput::DEFAULT_RSEM_PREP_REF = PATHS(RSEM_DEFAULT_EXE_DIR, RSEM_PREP_REF_EXE);
const std::string UserInput::DEFAULT_RSEM_CALC_EXP = PATHS(RSEM_DEFAULT_EXE_DIR, RSEM_CALC_EXP_EXE);
const std::string UserInput::DEFAULT_RSEM_CONV_SAM = PATHS(RSEM_DEFAULT_EXE_DIR, RSEM_CONV_SAM);

const std::string UserInput::GENEMARK_DEFAULT_EXE        = PATHS(FileSystem::get_exe_dir(),"/libs/gmst_linux_64/gmst.pl");
const std::string UserInput::DEFAULT_OUT_DIR             = PATHS(FileSystem::get_cur_dir(),"entap_outfiles");
const std::string UserInput::DEFAULT_INI_PATH            = PATHS(FileSystem::get_cur_dir(), ENTAP_INI_FILENAME);
const std::string UserInput::DEFAULT_ENTAP_DB_BIN_INI    = PATHS(FileSystem::get_exe_dir(), ENTAP_DATABASE_BIN_DEFAULT);
const std::string UserInput::DEFAULT_ENTAP_DB_SQL_INI    = PATHS(FileSystem::get_exe_dir(), ENTAP_DATABASE_SQL_DEFAULT);
const std::string UserInput::DEFAULT_EGG_SQL_DB_INI      = PATHS(FileSystem::get_exe_dir(), EGG_SQL_DB_DEFAULT);
const std::string UserInput::DEFAULT_EGG_DMND_DB_INI     = PATHS(FileSystem::get_exe_dir(), EGG_DMND_DEFAULT);
const std::string UserInput::DEFAULT_ENTAP_GRAPH_INI     = PATHS(FileSystem::get_exe_dir(), GRAPH_SCRIPT_DEF);
const std::string UserInput::DEFAULT_BUSCO_EXE           = "run_BUSCO.py";

// WARNING must match ENTAP_INPUT_FLAGS enum in UserInput.h
// If no default value exists, must use boost::any()
UserInput::EntapINIEntry UserInput::mUserInputs[] = {
        // Category      Input Flag                Shortened Input Flag   Description               Example         Variable Type           Default Value           Input Type          Parsed value
        {ENTAP_INI_NULL,ENTAP_INI_NULL           ,ENTAP_INI_NULL  ,ENTAP_INI_NULL             , ENTAP_INI_NULL  ,ENT_INI_VAR_STRING      ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},

/* General Input Commands */
        {INI_GENERAL   ,CMD_OUTPUT_DIR           ,ENTAP_INI_NULL  ,DESC_OUTPUT_DIR            ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,DEFAULT_OUT_DIR        ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_CONFIG               ,ENTAP_INI_NULL  ,DESC_CONFIG                ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_RUN_PROTEIN          ,ENTAP_INI_NULL  ,DESC_RUN_PROTEIN           ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_RUN_NUCLEO           ,ENTAP_INI_NULL  ,DESC_RUN_NUCLEO            ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_OVERWRITE            ,ENTAP_INI_NULL  ,DESC_OVERWRITE             ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_INI_FILE             ,ENTAP_INI_NULL  ,DESC_INI_FILE              ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,DEFAULT_INI_PATH       ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
//        {INI_GENERAL   ,CMD_HELP                 ,CMD_SHORT_HELP       ,DESC_HELP             ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL   ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
//        {INI_GENERAL   ,CMD_VERSION              ,CMD_SHORT_VERSION    ,DESC_VERSION          ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL   ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_INPUT_TRAN           ,CMD_SHORT_INPUT_TRAN ,DESC_INPUT_TRAN       ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_DATABASE             ,CMD_SHORT_DATABASE   ,DESC_DATABASE         ,ENTAP_INI_NULL   ,ENT_INI_VAR_MULTI_STRING,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_GRAPHING             ,ENTAP_INI_NULL  ,DESC_GRAPHING              ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_NO_TRIM              ,ENTAP_INI_NULL  ,DESC_NO_TRIM                  ,EX_NO_TRIM          ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_THREADS              ,CMD_SHORT_THREADS    ,DESC_THREADS          ,ENTAP_INI_NULL   ,ENT_INI_VAR_INT         ,DEFAULT_THREADS        ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_STATE                ,ENTAP_INI_NULL  ,DESC_STATE                 ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,DEFAULT_STATE          ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_NOCHECK              ,ENTAP_INI_NULL  ,DESC_NOCHECK               ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_GENERAL   ,CMD_OUTPUT_FORMAT        ,ENTAP_INI_NULL  ,DESC_OUTPUT_FORMAT         ,ENTAP_INI_NULL   ,ENT_INI_VAR_MULTI_INT   ,DEFAULT_OUT_FORMAT     ,ENT_INI_FILE         ,ENTAP_INI_NULL_VAL},

/* EnTAP Commands */
        {INI_ENTAP     ,CMD_ENTAP_DB_BIN         ,ENTAP_INI_NULL  ,DESC_ENTAP_DB_BIN          ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,DEFAULT_ENTAP_DB_BIN_INI, ENT_INI_FILE        ,ENTAP_INI_NULL_VAL},
        {INI_ENTAP     ,CMD_ENTAP_DB_SQL         ,ENTAP_INI_NULL  ,DESC_ENTAP_DB_SQL          ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,DEFAULT_ENTAP_DB_SQL_INI, ENT_INI_FILE        ,ENTAP_INI_NULL_VAL},
        {INI_ENTAP     ,CMD_ENTAP_GRAPH_PATH     ,ENTAP_INI_NULL  ,DESC_ENTAP_GRAPH_PATH      ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,DEFAULT_ENTAP_GRAPH_INI , ENT_INI_FILE        ,ENTAP_INI_NULL_VAL},
        {INI_ENTAP     ,ENTAP_INI_NULL           ,ENTAP_INI_NULL  ,ENTAP_INI_NULL             ,ENTAP_INI_NULL   ,ENT_INI_VAR_MULTI_INT,ENTAP_INI_NULL_VAL      , ENT_INPUT_FUTURE    ,ENTAP_INI_NULL_VAL},

/* Configuration Commands */
        {INI_CONFIG    ,CMD_DATA_GENERATE        ,ENTAP_INI_NULL  ,DESC_DATA_GENERATE         ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_CONFIG    ,CMD_DATABASE_TYPE        ,ENTAP_INI_NULL  ,DESC_DATABASE_TYPE         ,ENTAP_INI_NULL   ,ENT_INI_VAR_MULTI_INT   ,DEFAULT_DATA_TYPE      ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},

/* Expression Analysis Commands */
        {INI_EXPRESSION,CMD_FPKM                 ,ENTAP_INI_NULL  ,DESC_FPKM                  ,ENTAP_INI_NULL   ,ENT_INI_VAR_FLOAT       ,RSEM_FPKM_DEFAULT      ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},
        {INI_EXPRESSION,CMD_ALIGN_FILE           ,CMD_SHORT_ALIGN_FILE ,DESC_ALIGN_FILE       ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL},
        {INI_EXPRESSION,CMD_SINGLE_END           ,ENTAP_INI_NULL  ,DESC_SINGLE_END            ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},

/* Expression Analysis - RSEM Commands */
        {INI_EXP_RSEM  ,CMD_RSEM_CALC_EXP        ,ENTAP_INI_NULL  ,DESC_RSEM_CALC_EXP         ,EX_RSEM_CALC_EXP ,ENT_INI_VAR_STRING      ,DEFAULT_RSEM_CALC_EXP  ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},
        {INI_EXP_RSEM  ,CMD_RSEM_SAM_VALID       ,ENTAP_INI_NULL  ,DESC_RSEM_SAM_VALID        ,EX_RSEM_SAM_VALID,ENT_INI_VAR_STRING      ,DEFAULT_RSEM_SAM_VALID ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},
        {INI_EXP_RSEM  ,CMD_RSEM_PREP_REF        ,ENTAP_INI_NULL  ,DESC_RSEM_PREP_REF         ,EX_RSEM_PREP_REF ,ENT_INI_VAR_STRING      ,DEFAULT_RSEM_PREP_REF  ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},
        {INI_EXP_RSEM  ,CMD_RSEM_CONV_SAM        ,ENTAP_INI_NULL  ,DESC_RSEM_CONV_SAM         ,EX_RSEM_CONV_SAM ,ENT_INI_VAR_STRING      ,DEFAULT_RSEM_CONV_SAM  ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},

/* Frame Selection Commands */
        {INI_FRAME     ,CMD_COMPLETE_PROT        ,ENTAP_INI_NULL  ,DESC_COMPLETE_PROT         ,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL        ,ENTAP_INI_NULL_VAL     ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},
        {INI_FRAME     ,CMD_FRAME_SELECTION_FLAG ,ENTAP_INI_NULL  ,DESC_FRAME_SELECTION_FLAG  ,ENTAP_INI_NULL   ,ENT_INI_VAR_INT         ,DEFAULT_FRAME_SELECTION,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},

/* Frame Selection - GeneMarkST Commands */
        {INI_FRAME_GENEMARK,CMD_GENEMARKST_EXE   ,ENTAP_INI_NULL  ,DESC_GENEMARKST_EXE        ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,GENEMARK_DEFAULT_EXE   ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},

/* Frame Selection - TransDecoder Commands */
        {INI_FRAME_TRANSDECODER,CMD_TRANS_LONG_EXE,ENTAP_INI_NULL ,DESC_TRANS_LONG_EXE        ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,TRANSDECODER_LONG_DEFAULT_EXE, ENT_INI_FILE   ,ENTAP_INI_NULL_VAL},
        {INI_FRAME_TRANSDECODER,CMD_TRANS_PREDICT_EXE,ENTAP_INI_NULL,DESC_TRANS_PREDICT_EXE   ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,TRANSDECODER_PREDICT_DEFAULT_EXE, ENT_INI_FILE,ENTAP_INI_NULL_VAL},
        {INI_FRAME_TRANSDECODER,CMD_TRANS_MIN_FLAG,ENTAP_INI_NULL ,DESC_TRANS_MIN_FLAG        ,ENTAP_INI_NULL   ,ENT_INI_VAR_INT         ,DEFAULT_TRANSDECODER_MIN_PROTEIN, ENT_INI_FILE,ENTAP_INI_NULL_VAL},
        {INI_FRAME_TRANSDECODER,CMD_TRANS_NO_REF_START,ENTAP_INI_NULL, DESC_TRANS_NO_REF_START,ENTAP_INI_NULL   ,ENT_INI_VAR_BOOL       ,ENTAP_INI_NULL_VAL     ,ENT_INI_FILE          ,ENTAP_INI_NULL_VAL},

/* Similarity Search Commands */
        {INI_SIM_SEARCH,CMD_DIAMOND_EXE          ,ENTAP_INI_NULL  ,DESC_DIAMOND_EXE           ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,DIAMOND_DEFAULT_EXE    ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
        {INI_SIM_SEARCH,CMD_TAXON                ,ENTAP_INI_NULL  ,DESC_TAXON                 ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,ENTAP_INI_NULL_VAL     ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
        {INI_SIM_SEARCH,CMD_QCOVERAGE            ,ENTAP_INI_NULL  ,DESC_QCOVERAGE             ,ENTAP_INI_NULL   ,ENT_INI_VAR_FLOAT       ,DEFAULT_QCOVERAGE      ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
        {INI_SIM_SEARCH,CMD_TCOVERAGE            ,ENTAP_INI_NULL  ,DESC_TCOVERAGE             ,ENTAP_INI_NULL   ,ENT_INI_VAR_FLOAT       ,DEFAULT_TCOVERAGE      ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
        {INI_SIM_SEARCH,CMD_CONTAMINANT          ,CMD_SHORT_CONTAMINANT,DESC_CONTAMINANT      ,ENTAP_INI_NULL   ,ENT_INI_VAR_MULTI_STRING,ENTAP_INI_NULL_VAL     ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
        {INI_SIM_SEARCH,CMD_EVAL                 ,CMD_SHORT_EVAL       ,DESC_EVAL             ,ENTAP_INI_NULL   ,ENT_INI_VAR_FLOAT       ,DEFAULT_E_VALUE        ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
        {INI_SIM_SEARCH,CMD_UNINFORMATIVE        ,ENTAP_INI_NULL  ,DESC_UNINFORMATIVE         ,ENTAP_INI_NULL   ,ENT_INI_VAR_MULTI_STRING,DEFAULT_UNINFORMATIVE  ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},

/* Ontology Commands */
        {INI_ONTOLOGY  ,CMD_ONTOLOGY_FLAG        ,ENTAP_INI_NULL  ,DESC_ONTOLOGY_FLAG         ,ENTAP_INI_NULL   ,ENT_INI_VAR_MULTI_INT   ,DEFAULT_ONTOLOGY       ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
        {INI_ONTOLOGY  ,CMD_GO_LEVELS            ,ENTAP_INI_NULL  ,DESC_ONT_LEVELS            ,ENTAP_INI_NULL   ,ENT_INI_VAR_MULTI_INT   ,DEFAULT_ONT_LEVELS     ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},

/* Ontology - EggNOG Commands */
        {INI_ONT_EGGNOG,CMD_EGGNOG_SQL           ,ENTAP_INI_NULL  ,DESC_EGGNOG_SQL            ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,DEFAULT_EGG_SQL_DB_INI ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
        {INI_ONT_EGGNOG,CMD_EGGNOG_DMND          ,ENTAP_INI_NULL  ,DESC_EGGNOG_DMND           ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,DEFAULT_EGG_DMND_DB_INI,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
/* Ontology - InterPro Commands */
        {INI_ONT_INTERPRO,CMD_INTERPRO_EXE       ,ENTAP_INI_NULL  ,DESC_INTERPRO_EXE          ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,INTERPRO_DEF_EXE       ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},
        {INI_ONT_INTERPRO,CMD_INTER_DATA         ,ENTAP_INI_NULL  ,DESC_INTER_DATA            ,EX_INTER_DATA    ,ENT_INI_VAR_MULTI_STRING,ENTAP_INI_NULL_VAL     ,ENT_INI_FILE, ENTAP_INI_NULL_VAL},

/* Ontology - BUSCO Commands */
// DISABLED FOR NOW
        {INI_TRANSC_BUSCO,CMD_BUSCO_EXE          ,ENTAP_INI_NULL  ,DESC_BUSCO_EXE             ,EX_BUSCO_EXE     ,ENT_INI_VAR_STRING      ,DEFAULT_BUSCO_EXE      ,ENT_INPUT_FUTURE, ENTAP_INI_NULL_VAL},
        {INI_TRANSC_BUSCO,CMD_BUSCO_DATABASE     ,ENTAP_INI_NULL  ,DESC_BUSCO_DATABASE        ,ENTAP_INI_NULL   ,ENT_INI_VAR_STRING      ,ENTAP_INI_NULL_VAL     ,ENT_INPUT_FUTURE, ENTAP_INI_NULL_VAL},
        {INI_TRANSC_BUSCO,CMD_BUSCO_EVAL         ,ENTAP_INI_NULL  ,DESC_BUSCO_EVAL            ,ENTAP_INI_NULL   ,ENT_INI_VAR_FLOAT       ,DEFAULT_BUSCO_E_VALUE  ,ENT_INPUT_FUTURE, ENTAP_INI_NULL_VAL},

/* END COMMANDS */
        {ENTAP_INI_NULL,ENTAP_INI_NULL           ,ENTAP_INI_NULL  ,ENTAP_INI_NULL             , ENTAP_INI_NULL  ,ENT_INI_VAR_STRING      ,ENTAP_INI_NULL_VAL     ,ENT_COMMAND_LINE      ,ENTAP_INI_NULL_VAL}
};

UserInput::UserInput(int argc, const char** argv, FileSystem *fileSystem) {
    std::string ini_file_path;
    std::string root_dir;

    FS_dprint("Spawn Object - UserInput");

    mpFileSystem = fileSystem;

    // Parse command line arguments
    parse_arguments_tclap(argc, argv);

    if (has_input(INPUT_FLAG_OUTPUT_DIR)) {
        root_dir = get_user_input<ent_input_str_t>(INPUT_FLAG_OUTPUT_DIR);
    } else {
        root_dir = FileSystem::get_cur_dir();
    }
    mpFileSystem->set_root_dir(root_dir);

    ini_file_path = get_user_input<ent_input_str_t>(INPUT_FLAG_INI_FILE);

    // Ensure user has input the INI file path, EXIT otherwise
    if (!mpFileSystem->file_exists(ini_file_path)) {
        // INI file is required for EnTAP execution, generate one in the CWD
        ini_file_path = PATHS(mpFileSystem->get_cur_dir(), ENTAP_INI_FILENAME);
        generate_ini_file(ini_file_path);
        throw ExceptionHandler("INI file was not found and is required for EnTAP execution, generated at: " + ini_file_path,
                               ERR_ENTAP_CONFIG_CREATE_SUCCESS);
    } else {
        // Ini file exists and file path is valid
        mIniFilePath = ini_file_path;
        parse_ini(mIniFilePath);
    }
    parse_future_inputs();
}

UserInput::~UserInput() {
    FS_dprint("Killing object - UserInput");
}

/**
 * ======================================================================
 * Function void parse_future_inputs()
 *
 * Description          - Manages parsing and verifying of future inputs
 *                        that have defaults right now until they are used
 *                        as an actual input
 *                      - Will simply add these to user input database
 *
 * Notes                - Entry
 *
 *
 * @return              - None
 * ======================================================================
 */
void UserInput::parse_future_inputs() {
    ENTAP_INPUT_FLAGS input_flag;

    for (uint16 i=INPUT_FLAG_UNUSED; i <INPUT_FLAG_MAX; i++) {
        if (mUserInputs[i].input_type != ENT_INPUT_FUTURE) continue;

        input_flag = static_cast<ENTAP_INPUT_FLAGS>(i);

        switch (input_flag) {

            case INPUT_FLAG_ENTAP_HEADERS: {
                std::vector<ENTAP_HEADERS> header_vect(ENTAP_HEADER_COUNT-1);
                uint32 n = ENTAP_HEADER_UNUSED;
                std::generate(header_vect.begin(), header_vect.end(),
                              [&] { return static_cast<ENTAP_HEADERS>(++n); });
                mUserInputs[i].parsed_value = header_vect;
                break;
            }

            default:
                break;
        }
    }
}

/**
 * ======================================================================
 * Function void parse_ini(std::string &ini_path)
 *
 * Description          - Manages parsing and verifying the EnTAP ini file
 *                        provided by the user
 *
 * Notes                - Entry
 *
 * @param ini_path      - Path to the ini file
 *
 * @return              - None
 * ======================================================================
 */
void UserInput::parse_ini(std::string &ini_path) {
    FS_dprint("Parsing ini file at: " + ini_path);

    EntapINIEntry                               *ini_entry;
    std::string                                 line;
    std::string                                 key;
    std::string                                 val;

    if (!mpFileSystem->file_exists(ini_path)){
        throw ExceptionHandler("Ini file not found at: " + ini_path + " + exiting...", ERR_ENTAP_INPUT_PARSE);
    }

    // Parse ini file
    std::ifstream in_file(ini_path);
    while (std::getline(in_file,line)) {
        if (line.front() == INI_FILE_COMMENT) continue; // Skip INI file comments

        std::istringstream in_line(line);
        if (std::getline(in_line,key,INI_FILE_ASSIGN)) {
            // Ensure this INI file key is correct and user hasn't changed it EXIT otherwise
            ini_entry = check_ini_key(key);
            if (ini_entry == nullptr) {
                throw ExceptionHandler("Incorrect format in config file at line: " + in_line.str(), ERR_ENTAP_CONFIG_PARSE);
            } else {
                // Key is valid, parse the value
                if (std::getline(in_line,val)) {
                    // If value is empty, use default
                    if (val.empty()) {
                        ini_entry->parsed_value = ini_entry->default_value;
                    } else {
                        // Value is NOT empty, update the final value
                        try {
                            switch (ini_entry->var_type) {

                                case ENT_INI_VAR_INT:
                                    ini_entry->parsed_value = (ent_input_uint_t)std::stoi(val);
                                    break;
                                case ENT_INI_VAR_FLOAT:
                                    ini_entry->parsed_value = (ent_input_fp_t)std::stof(val);
                                    break;
                                case ENT_INI_VAR_STRING:
                                    ini_entry->parsed_value = val;
                                    break;
                                case ENT_INI_VAR_BOOL:
                                    if (val == INI_FILE_BOOL_TRUE) {
                                        ini_entry->parsed_value = true;
                                    } else if (val == INI_FILE_BOOL_FALSE){
                                        ;
                                    } else {
                                        throw ExceptionHandler("INI file boolean input must be true or false at line: "+ line,
                                            ERR_ENTAP_INPUT_PARSE);
                                    }
                                    break;
                                case ENT_INI_VAR_MULTI_FLOAT: {
                                    ent_input_multi_fp_t vect;
                                    std::stringstream ss(val);
                                    std::string token;
                                    while (std::getline(ss, token, INI_FILE_MULTI_DELIM)) {
                                        trim(token);
                                        vect.push_back(std::stof(token));
                                    }
                                    ini_entry->parsed_value = vect;
                                    break;
                                }
                                case ENT_INI_VAR_MULTI_INT: {
                                    ent_input_multi_int_t vect;
                                    std::stringstream ss(val);
                                    std::string token;
                                    while (std::getline(ss, token, INI_FILE_MULTI_DELIM)) {
                                        trim(token);
                                        vect.push_back((ent_input_uint_t)std::stoi(token));
                                    }
                                    ini_entry->parsed_value = vect;
                                    break;
                                }
                                case ENT_INI_VAR_MULTI_STRING: {
                                    ent_input_multi_str_t vect;
                                    std::stringstream ss(val);
                                    std::string token;
                                    while (std::getline(ss, token, INI_FILE_MULTI_DELIM)) {
                                        trim(token);
                                        vect.push_back(token);
                                    }
                                    ini_entry->parsed_value = vect;
                                    break;
                                }

                                default:
                                    break;
                            }
                        } catch (...) {
                            throw ExceptionHandler("Invalid INI file format at line: " + line, ERR_ENTAP_INPUT_PARSE);
                        }
                    }
                }
            }
        }
    }
    FS_dprint("Success!");
}

void UserInput::generate_ini_file(std::string &ini_path) {
    std::map<std::string, std::vector<EntapINIEntry*>> printed_categories;
    std::string line;

    std::ofstream ini_file(ini_path, std::ios::out | std::ios::trunc);

    for (EntapINIEntry& entry : mUserInputs) {
        // Skip NULL or future features
        if ((entry.category == ENTAP_INI_NULL) || (entry.input_type == ENT_INPUT_FUTURE)) continue;
        if (printed_categories.find(entry.category) != printed_categories.end()) {
            printed_categories[entry.category].push_back(&entry);
        } else {
            printed_categories.emplace(entry.category, std::vector<EntapINIEntry*>{&entry});
        }
    }

    try {

        // Print instructions on how to use the ini file
        ini_file << INI_FILE_COMMENT << "-------------------------------" << std::endl;
        ini_file << INI_FILE_COMMENT << " [" << "ini_instructions]"       << std::endl;
        ini_file << INI_FILE_COMMENT <<
                 "When using this ini file keep the following in mind:\n"              << INI_FILE_COMMENT <<
                  "\t1. Do not edit the input keys to the left side of the '=' sign\n" << INI_FILE_COMMENT <<
                  "\t2. Be sure to use the proper value type (either a string, list, or number)\n" << INI_FILE_COMMENT <<
                  "\t3. Do not add unecessary spaces to your input\n"                                         << INI_FILE_COMMENT <<
                  "\t4. When inputting a list, only add a '" << INI_FILE_MULTI_DELIM << "' between each entry" << std::endl;

        for (auto& pair : printed_categories) {

            // If this category has no entries, skip
            if (pair.second.empty()) {
                continue;
            }

            ini_file << INI_FILE_COMMENT << "-------------------------------" << std::endl;
            ini_file << INI_FILE_COMMENT << " [" << pair.first << "]" << std::endl;
            ini_file << INI_FILE_COMMENT << "-------------------------------" << std::endl;

            for (EntapINIEntry *entry : pair.second) {

                if (entry->input_type != ENT_INI_FILE) continue;

                // Print description
                std::stringstream ss(entry->description);
                while (std::getline(ss,line)) {
                    ini_file << INI_FILE_COMMENT << line << std::endl;
                }
                ss.str(""); ss.clear();

                // Print example
                ss.str(entry->example);
                while (std::getline(ss,line)) {
                    ini_file << INI_FILE_COMMENT << line << std::endl;
                }

                // Print input type
                ini_file << INI_FILE_COMMENT << "type:" << VAR_TYPE_STR[entry->var_type] << std::endl;

                ini_file << entry->input << INI_FILE_ASSIGN;

                switch (entry->var_type) {

                    case ENT_INI_VAR_BOOL:
                        if (entry->default_value.empty()) {
                            ini_file << INI_FILE_BOOL_FALSE;
                        } else {
                            ini_file << INI_FILE_BOOL_TRUE;
                        }
                        break;

                    case ENT_INI_VAR_STRING:
                        if (!entry->default_value.empty()) {
                            ini_file << boost::any_cast<ent_input_str_t>(entry->default_value);
                        }
                        break;
                    case ENT_INI_VAR_FLOAT:
                        if (!entry->default_value.empty()) {
                            ini_file << boost::any_cast<ent_input_fp_t >(entry->default_value);
                        }
                        break;
                    case ENT_INI_VAR_INT:
                        if (!entry->default_value.empty()) {
                            ini_file << boost::any_cast<ent_input_uint_t >(entry->default_value);
                        }
                        break;

                    case ENT_INI_VAR_MULTI_STRING: {
                        if (!entry->default_value.empty()) {
                            ent_input_multi_str_t vect = boost::any_cast<ent_input_multi_str_t>(entry->default_value);
                            if (vect.empty()) {
                                break;
                            } else {
                                for (auto &val : boost::any_cast<ent_input_multi_str_t>(entry->default_value)) {
                                    ini_file << val << INI_FILE_MULTI_DELIM;
                                }
                            }
                        }
                        break;
                    }
                    case ENT_INI_VAR_MULTI_INT: {
                        if (!entry->default_value.empty()) {
                            ent_input_multi_int_t vect = boost::any_cast<ent_input_multi_int_t>(entry->default_value);
                            if (vect.empty()) {
                                break;
                            } else {
                                for (auto &val : boost::any_cast<ent_input_multi_int_t>(entry->default_value)) {
                                    ini_file << std::to_string(val) << INI_FILE_MULTI_DELIM;
                                }
                            }
                        }
                        break;
                    }
                    case ENT_INI_VAR_MULTI_FLOAT: {
                        if (!entry->default_value.empty()) {
                            ent_input_multi_fp_t vect = boost::any_cast<ent_input_multi_fp_t>(entry->default_value);
                            if (vect.empty()) {
                                break;
                            } else {
                                for (auto &val : boost::any_cast<ent_input_multi_fp_t>(entry->default_value)) {
                                    ini_file << std::to_string(val) << INI_FILE_MULTI_DELIM;
                                }
                            }
                        }
                        break;
                    }

                    default:
                        continue;   // Continue if unrecognized VAR type
                }
                ini_file << std::endl;
            }
        }
    } catch (const std::exception &err) {
        ini_file.close();
        throw ExceptionHandler("Unable to generate INI file: " + std::string(err.what()),
                               ERR_ENTAP_CONFIG_CREATE);
    }
    ini_file.close();
}

/**
 * ======================================================================
 * Function bool check_ini_key(std::string&     key)
 *
 * Description          - Ensures EnTAP ini file has valid
 *                        entries and has not been edited by user
 *
 * Notes                - None
 *
 * @param key           - Key from configuration file
 * @return              - Flag if key is valid or not
 * =====================================================================
 */
UserInput::EntapINIEntry* UserInput::check_ini_key(std::string &key) {
    EntapINIEntry *ret = nullptr;

    LOWERCASE(key);
    for (EntapINIEntry &entry : mUserInputs) {
        if ((entry.input == key) && (entry.input_type == ENT_INI_FILE)) {
            ret = &entry;
            break;
        }
    }
    return ret;
}

/**
 * ======================================================================
 * Function void parse_arguments_tclap
 *                              (int            argc,
 *                               const char**   argv)
 *
 * Description          - Utilizes TCLAP libraries to parse user input
 *                        arguments
 *
 * Notes                - None
 *
 * @param argc          - User input size
 * @param argv          - User input
 *
 * @return              - None
 * ======================================================================
 */

 void UserInput::parse_arguments_tclap(int argc, const char ** argv) {
    std::vector<TCLAP::Arg*> tclap_arguments(INPUT_FLAG_MAX);
    TCLAP::Arg *arg = nullptr;
    EntapINIEntry *entry = nullptr;
    boost::any any_val;

    try {
        TCLAP::CmdLine cmd("EnTAP\nAlexander Hart and Dr. Jill Wegrzyn\nUniversity of Connecticut\nCopyright 2017-2021",
                           ' ', ENTAP_VERSION_STR);

        // Generate Arguments and add them to CMD
        for (uint32 i = 1; i < INPUT_FLAG_MAX; i++) {

            tclap_arguments[i] = nullptr;
            entry = &mUserInputs[i];
            // Continue if entry is NOT parsed through command line
            if (entry->input_type != ENT_COMMAND_LINE) continue;    // CONTINUE

            if (!entry->default_value.empty()) {
                any_val = entry->default_value;
            }

            switch (entry->var_type) {

                // These are considered Switch Arguments
                case ENT_INI_VAR_BOOL:
                    arg = new TCLAP::SwitchArg(entry->short_input,entry->input, entry->description, cmd, false);
                    break;

                // These are considered Value Arguments
                case ENT_INI_VAR_STRING:
                    if (entry->default_value.empty()) {
                        any_val = ENTAP_INI_NULL;
                    }
                    arg = new TCLAP::ValueArg<ent_input_str_t>(entry->short_input, entry->input, entry->description, false,
                                                          boost::any_cast<ent_input_str_t>(any_val), "string",cmd);
                    break;

                case ENT_INI_VAR_INT:
                    if (entry->default_value.empty()) {
                        any_val = ENTAP_INI_NULL_INT;
                    }
                    arg = new TCLAP::ValueArg<ent_input_uint_t >(entry->short_input, entry->input, entry->description, false,
                                                          boost::any_cast<ent_input_uint_t>(any_val), "integer", cmd);
                    break;

                case ENT_INI_VAR_FLOAT:
                    if (entry->default_value.empty()) {
                        any_val = ENTAP_INI_NULL_FLT;
                    }
                    arg = new TCLAP::ValueArg<ent_input_fp_t >(entry->short_input, entry->input, entry->description, false,
                                                     boost::any_cast<ent_input_fp_t>(any_val), "decimal", cmd);
                    break;

                // These are considered Multi Arguments
                case ENT_INI_VAR_MULTI_STRING:
                    arg = new TCLAP::MultiArg<ent_input_str_t>(entry->short_input, entry->input, entry->description, false, "string list", cmd);
                    break;

                case ENT_INI_VAR_MULTI_INT:
                    arg = new TCLAP::MultiArg<ent_input_uint_t> (entry->short_input, entry->input, entry->description, false, "integer list", cmd);
                    break;

                case ENT_INI_VAR_MULTI_FLOAT:
                    arg = new TCLAP::MultiArg<ent_input_fp_t>(entry->short_input, entry->input, entry->description, false, "decimal list", cmd);
                    break;

                // Skip unhandled types
                default:
                    break;
            }
            tclap_arguments[i] = arg;
        }

        // Parse the command line
        cmd.parse( argc, argv );

        // if no inputs, print usage and exit
        if (argc == 1) {
            TCLAP::StdOutput out;
            out.usage(cmd);

            for (auto &argument : tclap_arguments) {
                delete argument;
            }

            throw ExceptionHandler("",0);   // EXIT NO ARGUMENTS
        }

        // Loop through arguments and populate our overall data
        for (uint32 i = 0; i < tclap_arguments.size(); i++) {
            arg = tclap_arguments[i];
            entry = &mUserInputs[i];

            // Ensure it is a command line argument
            if (entry->input_type == ENT_COMMAND_LINE) {

                // If argument NULL or not set by the user
                if (arg == nullptr || !arg->isSet()) {
                    // Yes, set the final value to the default value
                    entry->parsed_value = entry->default_value;
                } else {
                    // No, we want to update final value with the user input value
                    switch (entry->var_type) {

                        // If switch argument
                        case ENT_INI_VAR_BOOL:
                            entry->parsed_value = true;
                            break;

                        // If value argument
                        case ENT_INI_VAR_STRING: {
                            auto *value_arg = dynamic_cast<TCLAP::ValueArg<ent_input_str_t>*>(arg);
                            entry->parsed_value = value_arg->getValue();
                            break;
                        }

                        case ENT_INI_VAR_INT: {
                            auto *value_arg = dynamic_cast<TCLAP::ValueArg<ent_input_uint_t > *>(arg);
                            entry->parsed_value = value_arg->getValue();
                            break;
                        }

                        case ENT_INI_VAR_FLOAT: {
                            auto *value_arg = dynamic_cast<TCLAP::ValueArg<ent_input_fp_t > *>(arg);
                            entry->parsed_value = value_arg->getValue();
                            break;
                        }

                        // If multi argument
                        case ENT_INI_VAR_MULTI_STRING: {
                            auto *multi_arg = dynamic_cast<TCLAP::MultiArg<ent_input_str_t> *>(arg);
                            entry->parsed_value = multi_arg->getValue();
                            break;
                        }

                        case ENT_INI_VAR_MULTI_INT: {
                            auto *multi_arg = dynamic_cast<TCLAP::MultiArg<ent_input_uint_t> *>(arg);
                            entry->parsed_value = multi_arg->getValue();
                            break;
                        }
                        case ENT_INI_VAR_MULTI_FLOAT: {
                            auto *multi_arg = dynamic_cast<TCLAP::MultiArg<ent_input_fp_t> *>(arg);
                            entry->parsed_value = multi_arg->getValue();
                            break;
                        }

                        default:
                            break;
                    }

                }

            } else {
                ;
            }

        }
        for (auto &argument : tclap_arguments) {
            delete argument;
        }

    } catch (TCLAP::ArgException &e) {
        for (auto &argument : tclap_arguments) {
            delete argument;
        }
        throw ExceptionHandler(e.what(), ERR_ENTAP_INPUT_PARSE);
    }
 }

/**
 * ======================================================================
 * Function void verify_user_input(void)
 *
 * Description          - Performs sanity checks on user input
 *
 * Notes                - None
 *
 * @return              - None
 * =====================================================================
 */
bool UserInput::verify_user_input() {

    bool                     is_interpro;
    bool                     is_protein;
    bool                     is_nucleotide;
    bool                     is_config;
    bool                     is_run;
    std::string              species;
    std::string              input_tran_path;
    std::vector<uint16>      ont_flags;
    ent_input_multi_int_t    GO_levels;
    EntapDatabase           *pEntap_database = nullptr;
    std::unique_ptr<QueryData> pQuery_Data=nullptr;

    // If graphing flag, check if it is allowed then EXIT
    if (has_input(INPUT_FLAG_GRAPH)) {
        std::string graphing_exe = get_user_input<ent_input_str_t>(INPUT_FLAG_ENTAP_GRAPH);
        if (!mpFileSystem->file_exists(graphing_exe)) {
            std::cout<<"Graphing is NOT enabled on this system! Graphing script could not "
                    "be found at: "<<graphing_exe << std::endl;
        }
        GraphingManager gmanager = GraphingManager(graphing_exe, mpFileSystem);
        if (gmanager.is_graphing_enabled()) {
            std::cout<< "Graphing is enabled on this system!" << std::endl;
            throw ExceptionHandler("",ERR_ENTAP_SUCCESS);
        } else {
            std::cout<<"Graphing is NOT enabled on this system!,"
                    " ensure that you have python with the Matplotlib module installed."<<std::endl;
            throw ExceptionHandler("",ERR_ENTAP_SUCCESS);
        }
    }


    // ------------ Config / Run Required beyond this point ---------------- //

    is_config     = has_input(INPUT_FLAG_CONFIG);     // ignore 'config config'
    is_protein    = has_input(INPUT_FLAG_RUNPROTEIN);
    is_nucleotide = has_input(INPUT_FLAG_RUNNUCLEOTIDE);

    mIsConfig = is_config;

    if (is_protein && is_nucleotide) {
        throw ExceptionHandler("Cannot specify both protein and nucleotide input flags",
                               ERR_ENTAP_INPUT_PARSE);
    }
    is_run = is_protein || is_nucleotide;

    // Check config and run flags
    if (!is_config && !is_run) {
        throw(ExceptionHandler("Either config option or run option are required",
                               ERR_ENTAP_INPUT_PARSE));
    } else if (is_config && is_run) {
        throw(ExceptionHandler("Cannot specify both config and run flags",
                               ERR_ENTAP_INPUT_PARSE));
    }
    print_user_input();

    // If user wants to skip this check, EXIT (do this after we print the input for debugging purposes)
    if (has_input(INPUT_FLAG_NOCHECK)) {
        FS_dprint("WARNING User is skipping input verification!! :(");
        return is_config;
    }

    try {

        verify_databases(is_run);

        // Handle generic flags
        if (has_input(INPUT_FLAG_OUTPUT_FORMAT)) {
            ent_input_multi_int_t output_formats = get_user_input<ent_input_multi_int_t>(INPUT_FLAG_OUTPUT_FORMAT);
            for (ent_input_uint_t flag : output_formats) {
                if (flag <= FileSystem::ENT_FILE_UNUSED || flag >= FileSystem::ENT_FILE_OUTPUT_FORMAT_MAX) {
                    throw ExceptionHandler("Invalid flag for Output Format (" + std::to_string(flag) + ")",
                                           ERR_ENTAP_INPUT_PARSE);
                }
            }
        }

        // Handle EnTAP execution commands
        if (is_run) {

            // Verify EnTAP database can be generated
            FS_dprint("Verifying EnTAP database...");
            pEntap_database = new EntapDatabase(mpFileSystem);
            // Find database type that will be used by the rest (use 0 index no matter what)
            ent_input_multi_int_t entap_database_types =
                    get_user_input<ent_input_multi_int_t>(INPUT_FLAG_DATABASE_TYPE);
            auto type = static_cast<EntapDatabase::DATABASE_TYPE>(entap_database_types[0]);
            if (!pEntap_database->set_database(type, get_entap_database_path(type))) {
                throw ExceptionHandler("Unable to open EnTAP database from paths given" + pEntap_database->print_error_log(),
                                       ERR_ENTAP_READ_ENTAP_DATA_GENERIC);
            }
            // Verify database type
            if (!pEntap_database->is_valid_version()) {
                throw ExceptionHandler("EnTAP database version incompatible with this version of software. Please download the correct version.\nYou have: " +
                                               pEntap_database->get_current_version_str() + "\nYou need: " +
                                               pEntap_database->get_required_version_str(), ERR_ENTAP_READ_ENTAP_DATA_GENERIC);
            }

            FS_dprint("Success!");

            // Verify input transcriptome
            if (!has_input(INPUT_FLAG_TRANSCRIPTOME)) {
                throw(ExceptionHandler("Must enter a valid transcriptome",ERR_ENTAP_INPUT_PARSE));
            } else {
                input_tran_path = get_user_input<ent_input_str_t>(INPUT_FLAG_TRANSCRIPTOME);
                if (!mpFileSystem->file_exists(input_tran_path)) {
                    throw(ExceptionHandler("Transcriptome not found at: " + input_tran_path,
                                           ERR_ENTAP_INPUT_PARSE));
                } else if (mpFileSystem->file_empty(input_tran_path)) {
                    throw(ExceptionHandler("Transcriptome file empty: "+ input_tran_path,
                                           ERR_ENTAP_INPUT_PARSE));
                } else if (!mpFileSystem->check_fasta(input_tran_path)) {
                    throw(ExceptionHandler("File not in fasta format or corrupt! "+ input_tran_path,
                                           ERR_ENTAP_INPUT_PARSE));
                } else {
                    // File valid, verify that the transcriptome is in correct format
                    try {
                        pQuery_Data = std::unique_ptr<QueryData>(new QueryData(input_tran_path, this, mpFileSystem));
                    } catch (const ExceptionHandler &e) {
                        throw e;
                    }

                }
            }

            // Verify species for taxonomic relevance
            if (has_input(INPUT_FLAG_SPECIES)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_SPECIES].input);
                verify_species(SPECIES, pEntap_database);
            }

            // Verify contaminant
            if (has_input(INPUT_FLAG_CONTAMINANT)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_CONTAMINANT].input);
                verify_species(CONTAMINANT, pEntap_database);
            }

            // Verify path + extension for alignment file
            if (has_input(INPUT_FLAG_ALIGN)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_ALIGN].input);
                ent_input_str_t align_file = get_user_input<ent_input_str_t>(INPUT_FLAG_ALIGN);
                std::string align_ext = mpFileSystem->get_file_extension(align_file, false);
                std::transform(align_ext.begin(), align_ext.end(), align_ext.begin(), ::tolower);
                if (!mpFileSystem->file_exists(align_file)) {
                    throw ExceptionHandler("BAM/SAM file not found at: " + align_file + " exiting...",
                                           ERR_ENTAP_INPUT_PARSE);
                }
                if (align_ext != FileSystem::EXT_SAM && align_ext != FileSystem::EXT_BAM) {
                    throw ExceptionHandler("Alignment file must have a .bam or .sam extension",
                                           ERR_ENTAP_INPUT_PARSE);
                }
            }

            // Verify FPKM
            if (has_input(INPUT_FLAG_FPKM)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_FPKM].input);
                ent_input_fp_t fpkm = get_user_input<ent_input_fp_t>(INPUT_FLAG_FPKM);
                if (fpkm > FPKM_MAX || fpkm < FPKM_MIN) {
                    throw ExceptionHandler("Selected FPKM threshold is out of range, must be between " + std::to_string(FPKM_MIN) +
                                           " and " + std::to_string(FPKM_MAX), ERR_ENTAP_INPUT_PARSE);
                }
            }

            // Verify query coverage
            if (has_input(INPUT_FLAG_QCOVERAGE)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_QCOVERAGE].input);
                ent_input_fp_t qcoverage = get_user_input<ent_input_fp_t>(INPUT_FLAG_QCOVERAGE);
                if (qcoverage > COVERAGE_MAX || qcoverage < COVERAGE_MIN) {
                    throw ExceptionHandler("Query coverage is out of range, but be between " +
                                           std::to_string(COVERAGE_MIN) +
                                           " and " + std::to_string(COVERAGE_MAX), ERR_ENTAP_INPUT_PARSE);
                }
            }

            // Verify target coverage
            if (has_input(INPUT_FLAG_TCOVERAGE)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_TCOVERAGE].input);
                ent_input_fp_t qcoverage = get_user_input<ent_input_fp_t>(INPUT_FLAG_TCOVERAGE);
                if (qcoverage > COVERAGE_MAX || qcoverage < COVERAGE_MIN) {
                    throw ExceptionHandler("Target coverage is out of range, but be between " +
                                           std::to_string(COVERAGE_MIN) +
                                           " and " + std::to_string(COVERAGE_MAX), ERR_ENTAP_INPUT_PARSE);
                }
            }

            // Verify Ontology Flags
            is_interpro = false;
            if (has_input(INPUT_FLAG_ONTOLOGY)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_ONTOLOGY].input);
                ont_flags = get_user_input<ent_input_multi_int_t>(INPUT_FLAG_ONTOLOGY);
                for (ent_input_uint_t i = 0; i < ont_flags.size() ; i++) {
                    if ((ont_flags[i] > ONT_SOFTWARE_COUNT) || ont_flags[i] < 0) {
                        throw ExceptionHandler("Invalid ontology flags being used", ERR_ENTAP_INPUT_PARSE);
                    }
                    if (ont_flags[i] == ONT_INTERPRO_SCAN && !is_interpro) is_interpro = true;
                }
            }

            // Verify Go Levels
            if (has_input(INPUT_FLAG_GO_LEVELS)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_GO_LEVELS].input);
                GO_levels = get_user_input<ent_input_multi_int_t>(INPUT_FLAG_GO_LEVELS);

                if (GO_levels.size() > MAX_GO_LEVELS_SELECTED) {
                    throw ExceptionHandler("Only up to " + std::to_string(MAX_GO_LEVELS_SELECTED) + " Gene Ontology"
                          " levels can be selected", ERR_ENTAP_INPUT_PARSE);
                }

                for (uint16 lvl : GO_levels) {
                    if (lvl >= MAX_GO_LEVEL) {
                        throw ExceptionHandler("The maximum Gene Ontology level that can be selected is " +
                            std::to_string(MAX_GO_LEVEL), ERR_ENTAP_INPUT_PARSE);
                    }
                }
            }

            // Verify InterPro databases
            if (is_interpro && !ModInterpro::valid_input(this)) {
                throw ExceptionHandler("InterProScan selected, but invalid databases input!", ERR_ENTAP_INPUT_PARSE);
            }

            // Verify Frame Selection flag
            if (has_input(INPUT_FLAG_FRAME_SELECTION)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_FRAME_SELECTION].input);
                ent_input_uint_t frame_software = get_user_input<ent_input_uint_t >(INPUT_FLAG_FRAME_SELECTION);
                if (frame_software <= FRAME_UNUSED || frame_software >= FRAME_SOFTWARE_COUNT) {
                    throw ExceptionHandler("Invalid Frame Selection software input!", ERR_ENTAP_INPUT_PARSE);
                }
            }

            // Verify minimum protein value
            if (has_input(INPUT_FLAG_TRANS_MIN_PROTEIN)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_TRANS_MIN_PROTEIN].input);
                ent_input_uint_t  min_protein = get_user_input<ent_input_uint_t>(INPUT_FLAG_TRANS_MIN_PROTEIN);
                if (min_protein < TRANS_MIN_PROTEIN_MIN) {
                    throw ExceptionHandler("Invalid TransDecoder minimum protein length value, minimum is: " +
                        std::to_string(TRANS_MIN_PROTEIN_MIN), ERR_ENTAP_INPUT_PARSE);
                }
            }

        } else {
            // Must be config

            // Verify BUSCO database if the user has input it
            if (has_input(INPUT_FLAG_BUSCO_DATABASE)) {
                FS_dprint("Verifying input flag " + mUserInputs[INPUT_FLAG_BUSCO_DATABASE].input);
                ent_input_multi_str_t busco_db = get_user_input<ent_input_multi_str_t>(INPUT_FLAG_BUSCO_DATABASE);
                BuscoDatabase buscoDatabase = BuscoDatabase(mpFileSystem);
                for (std::string &database : busco_db) {
                    std::string temp;
                    if (!buscoDatabase.valid_database(database, temp)) {
                        throw ExceptionHandler("Invalid BUSCO database entered must be a URL or valid database name: " + database,
                                               ERR_ENTAP_INPUT_PARSE);
                    }
                }

            }
        }

        // Verify software path for both CONFIG and EXECUTION
        if (get_user_input<ent_input_str_t>(INPUT_FLAG_STATE) == DEFAULT_STATE) {
            std::string state = DEFAULT_STATE;
            // only handling default now
            verify_software_paths(state, is_protein, is_run,pQuery_Data.get());
        }

    }catch (const ExceptionHandler &e) {
        SAFE_DELETE(pEntap_database);
        throw e;
    }
    FS_dprint("Success! Input verified");
    SAFE_DELETE(pEntap_database);
    return is_config;
}


/**
 * ======================================================================
 * Function void verify_databases(bool is_run)
 *
 * Description          - Ensures the user is entering valid databases
 *                        (paths exist, DIAMOND extension if needed)
 *
 * Notes                - None
 *
 * @param is_run        - True is we are running, false if configuration
 * @return              - None
 * ======================================================================
 */
void UserInput::verify_databases(bool is_run) {

    ent_input_multi_str_t     other_data;

    if (has_input(INPUT_FLAG_DATABASE)) {
        other_data = get_user_input<ent_input_multi_str_t>(INPUT_FLAG_DATABASE);
    } else if (is_run){
        // Must specify database when executing
        throw ExceptionHandler("Must select databases when executing main pipeline", ERR_ENTAP_INPUT_PARSE);
    }

    if (other_data.size() > MAX_DATABASE_SIZE) {
        throw ExceptionHandler("Too many databases selected, the max is " + std::to_string(MAX_DATABASE_SIZE), ERR_ENTAP_INPUT_PARSE);
    }

    // Check each database entered
    for (auto const& path: other_data) {
        if (!mpFileSystem->file_exists(path) || mpFileSystem->file_empty(path)) {
            throw ExceptionHandler("Database path invalid or empty: " + path, ERR_ENTAP_INPUT_PARSE);
        }
        FS_dprint("User has input a database at: " + path);
        // Is file extension diamond?
        if (mpFileSystem->get_file_extension(path, false) == FileSystem::EXT_DMND) {
            // Yes, are we configuring?
            if (!is_run) {
                throw ExceptionHandler("Cannot input DIAMOND (.dmnd) database when configuring!", ERR_ENTAP_INPUT_PARSE);
            }
        } else {
            // No, are we executing main pipeline?
            if (is_run) {
                throw ExceptionHandler("Must input DIAMOND (.dmnd) database when executing!", ERR_ENTAP_INPUT_PARSE);
            }
        }
    }
}


/**
 * ======================================================================
 * Function void print_user_input()
 *
 * Description          - Handles printing of user selected flags to
 *                        EnTAP statistics/log file
 *                      - All execution paths are NOT required
 *
 * Notes                - Accesses global software execution paths from
 *                        config file
 *
 * @return              - None
 *
 * =====================================================================
 */
void UserInput::print_user_input() {


    std::string         output;
    std::stringstream   ss;
    std::string         config_text;
    std::string         key;
    ENTAP_INPUT_FLAGS   input_flag;
    std::time_t         time;
    std::chrono::time_point<std::chrono::system_clock> start_time;

    FS_dprint("Printing user input...");
    start_time = std::chrono::system_clock::now();
    time = std::chrono::system_clock::to_time_t(start_time);
    mIsConfig ? config_text = "Configuration" : config_text = "Execution";

    mpFileSystem->format_stat_stream(ss, "EnTAP Run Information - " + config_text);

    ss <<
       "Current EnTAP Version: "   << ENTAP_VERSION_STR            <<
       "\nStart time: "            << std::ctime(&time)            <<
       "\nWorking directory has been set to: "  << mpFileSystem->get_root_path()<<
       "\nUser Inputs:";

    // Print all user inputs (fairly raw right now)

    for (uint32 i=1; i < INPUT_FLAG_MAX; i++) {

        const auto* it = &mUserInputs[i];

        if (it == nullptr || (it->input_type == ENT_INPUT_FUTURE)) continue;

        key = it->input;
        input_flag = static_cast<ENTAP_INPUT_FLAGS>(i);
        ss << "\n" << key << ": ";

        // Skip entering a value if it is NULL or BOOL (we want to print true or false here)
        if (!has_input(input_flag) && it->var_type != ENT_INI_VAR_BOOL) {
            continue;
        }

        switch (it->var_type) {

            case ENT_INI_VAR_BOOL:
                if (!has_input(input_flag)) {
                    ss << INI_FILE_BOOL_FALSE;
                } else {
                    ss << INI_FILE_BOOL_TRUE;
                }
                break;

            case ENT_INI_VAR_STRING:
                ss << get_user_input<ent_input_str_t>(input_flag);
                break;

            case ENT_INI_VAR_INT: {
                auto v = get_user_input<ent_input_uint_t>(input_flag);
                ss << std::to_string(v);
                break;
            }

            case ENT_INI_VAR_FLOAT: {
                auto v = get_user_input<ent_input_fp_t>(input_flag);
                ss << std::to_string(v);
                break;
            }

            case ENT_INI_VAR_MULTI_STRING: {
                auto v = get_user_input<ent_input_multi_str_t>(input_flag);
                for (auto &ii : v) {
                    ss << ii << INI_FILE_MULTI_DELIM;
                }
                break;
            }

            case ENT_INI_VAR_MULTI_INT: {
                auto v = get_user_input<ent_input_multi_int_t>(input_flag);
                for (auto &ii : v) {
                    ss << std::to_string(ii) << INI_FILE_MULTI_DELIM;
                }
                break;
            }

            case ENT_INI_VAR_MULTI_FLOAT: {
                auto v = get_user_input<ent_input_multi_fp_t>(input_flag);
                for (auto &ii : v) {
                    ss << std::to_string(ii) << INI_FILE_MULTI_DELIM;
                }
                break;
            }

            default:
                continue;
        }
    }
    // Print to log file and debug
    output = ss.str() + "\n";
    mpFileSystem->print_stats(output);
    FS_dprint(output+"\n");
}


/**
 * ======================================================================
 * Function void verify_species(SPECIES_FLAGS flag, EntapDatabase *database)
 *
 * Description          - Verify species/tax level input by the user
 *                      - Ensure it can be found within the tax database
 *                      - Ensure it's in the right format
 *
 * Notes                - None
 *
 * @param exe           - Boost map of user inputs
 * @return              - None
 * ======================================================================
 */
void UserInput::verify_species(SPECIES_FLAGS flag, EntapDatabase *database) {

    ent_input_multi_str_t    species;
    std::string              raw_species;

    if (database == nullptr) return;

    if (flag == SPECIES) {
        raw_species = get_target_species_str();
        species.push_back(raw_species);
    } else if (flag == CONTAMINANT) {
        species = get_user_input<ent_input_multi_str_t>(INPUT_FLAG_CONTAMINANT);
        for (std::string &contam : species) {
            process_user_species(contam);
        }
    }
    if (species.empty()) return;

    for (std::string &s : species) {
        if (database->get_tax_entry(s).is_empty()) {
            throw ExceptionHandler("Error in one of your inputted taxons: " + s + " it is not located"
                                   " within the taxonomic database. You may remove it or select another",
                                    ERR_ENTAP_INPUT_PARSE);
        } else {
            FS_dprint("Verified species: " + s);
        }
    }
}

/**
 * ======================================================================
 * Function void process_user_species(std::string &input)
 *
 * Description          - Format species user has input
 *
 * Notes                - Throw error on failure
 *
 * @param input         - Species to be formatted
 * @return              - None
 * ======================================================================
 */
void UserInput::process_user_species(std::string &input) {
    LOWERCASE(input);
    STR_REPLACE(input, '_', ' ');
}

/**
 * ======================================================================
 * Function void verify_state(std::string &state, bool runP, bool is_execution
 *                            std::vector<uint16> &ontology)
 *
 * Description          - Entry to check execution paths for software based
 *                        on state
 *                      - Sanity check on software specific commands
 *
 * Notes                - Throw error on failure
 *
 * @param state         - State inputted by user (or default)
 * @param runP          - Blastp flag (yes/no)
 * @param is_execution  - If we are running execution or configuration
 * @param ontology      - Vector of ontology flags
 *
 * @return              - None
 * ======================================================================
 */
void UserInput::verify_software_paths(std::string &state, bool runP, bool is_execution, QueryData *pQuery_data) {
    uint8 execute = 0x0;
    std::pair<bool, std::string> out;
    ent_input_str_t dmnd_exe;
    ent_input_str_t egg_db_dmnd;
    ent_input_str_t egg_db_sql;
    ent_input_str_t interpro_exe;
    ent_input_str_t busco_exe;
    ent_input_str_t genemarkst_exe;
    ent_input_str_t transdecoder_exe_predict;
    ent_input_str_t transdecoder_exe_longorf;
    ent_input_multi_int_t ontology_flags;
    ent_input_uint_t frame_selection_software;

    dmnd_exe     = get_user_input<ent_input_str_t>(INPUT_FLAG_DIAMOND_EXE);
    egg_db_dmnd  = get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_DMND_DB);
    egg_db_sql   = get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_SQL_DB);
    interpro_exe = get_user_input<ent_input_str_t>(INPUT_FLAG_INTERPRO_EXE);
    busco_exe    = get_user_input<ent_input_str_t>(INPUT_FLAG_BUSCO_EXE);
    genemarkst_exe = get_user_input<ent_input_str_t>(INPUT_FLAG_GENEMARKST_EXE);
    transdecoder_exe_predict = get_user_input<ent_input_str_t>(INPUT_FLAG_TRANS_PREDICT_EXE);
    transdecoder_exe_longorf = get_user_input<ent_input_str_t>(INPUT_FLAG_TRANS_LONGORF_EXE);

    ontology_flags = get_user_input<ent_input_multi_int_t>(INPUT_FLAG_ONTOLOGY);
    frame_selection_software = get_user_input<ent_input_uint_t>(INPUT_FLAG_FRAME_SELECTION);

    if (state == DEFAULT_STATE) {
        bool is_run_frame_select;

        if (run_expression_filtering()) {
            execute |= EXPRESSION_FILTERING;
        }

        if (run_frame_selection(pQuery_data, is_run_frame_select)) {
            if (is_run_frame_select) {
                execute |= FRAME_SELECTION;
            }
        }
        execute |= SIMILARITY_SEARCH;
        execute |= GENE_ONTOLOGY;

    }
    FS_dprint("Verifying software...");

    if (is_execution) {

        // add expression filtering

        // Check frame selection software
        if (execute & FRAME_SELECTION) {
            switch (frame_selection_software) {

                case FRAME_GENEMARK_ST:
                    FS_dprint("Verifying that GeneMarkS-T is executable...");
                    if (!ModGeneMarkST::is_executable(genemarkst_exe)) {
                        throw ExceptionHandler("Could not execute a test run of GeneMarkS-T, be sure "
                                               "it's properly installed and the executable is correct",
                                               ERR_ENTAP_INPUT_PARSE);
                    }
                    break;

                case FRAME_TRANSDECODER:
                    FS_dprint("Verifying that Transdecoder is executable");
                    if (!ModTransdecoder::is_executable(transdecoder_exe_longorf, transdecoder_exe_predict)) {
                        throw ExceptionHandler("Could not execute a test run of Transdecoder, be sure "
                                               "it's properly installed and the executable is correct",
                                               ERR_ENTAP_INPUT_PARSE);
                    }
                    break;
            }
        }

        // Check SIMLARITY SEARCH software
        if (execute & SIMILARITY_SEARCH) {
            if (!ModDiamond::is_executable(dmnd_exe)) {
                throw ExceptionHandler("Could not execute a test run of DIAMOND, be sure it's properly "
                                       "installed and the path is correct", ERR_ENTAP_INPUT_PARSE);
            }
        }
        if (execute & GENE_ONTOLOGY) {
            for (uint16 flag : ontology_flags) {
                switch (flag) {
#ifdef EGGNOG_MAPPER
                    case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
                    if (!pFileSystem->file_exists(EGG_SQL_DB_PATH))
                        return std::make_pair(false, "Could not find EggNOG SQL database at: " + EGG_SQL_DB_PATH);
                    if (!ModEggnog::is_executable())
                        return std::make_pair(false, "Test of EggNOG Emapper failed, "
                                "ensure python is properly installed and the paths are correct");
                    break;
#endif
                    case ONT_INTERPRO_SCAN:
                        FS_dprint("Verifying InterProScan inputs...");
                        if (!ModInterpro::is_executable(interpro_exe)) {
                            throw ExceptionHandler("Could not execute test run of InterProScan with execution command: " +
                                interpro_exe, ERR_ENTAP_INPUT_PARSE);
                        }
                        FS_dprint("Success!");
                        break;

                    case ONT_EGGNOG_DMND:
                        FS_dprint("Verifying EggNOG inputs...");
                        if (!mpFileSystem->file_exists(egg_db_sql))
                            throw ExceptionHandler("Could not find EggNOG SQL database at: " + egg_db_sql, ERR_ENTAP_INPUT_PARSE);
                        else if (!mpFileSystem->file_exists(egg_db_dmnd))
                            throw ExceptionHandler("Could not find EggNOG Diamond Database at: " + egg_db_dmnd, ERR_ENTAP_INPUT_PARSE);
                        else if (!ModEggnogDMND::is_executable(dmnd_exe))
                            throw ExceptionHandler("Could not execute a test run of DIAMOND for EggNOG analysis has failed", ERR_ENTAP_INPUT_PARSE);
                        FS_dprint("Success!");
                        break;

                    case ONT_BUSCO:
                        FS_dprint("Verifying BUSCO inputs...");
                        if (ModBUSCO::is_executable(busco_exe)) {
                            ;
                        } else {
                            throw ExceptionHandler("Could not execute test run of BUSCO with execution command: " +
                                busco_exe, ERR_ENTAP_INPUT_PARSE);
                        }
                        FS_dprint("Success!");
                        break;

                    default:
                        break;
                }
            }
        }
    // No, using CONFIGURATION stage of pipeline
    } else {

        // Check if EggNOG DIAMOND database exists, if not, check DIAMOND run
        if (!mpFileSystem->file_exists(get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_DMND_DB))) {
            if (!ModDiamond::is_executable(dmnd_exe)) {
                throw ExceptionHandler("EggNOG DIAMOND database was not found at: " + egg_db_dmnd +
                                       "\nThe DIAMOND test run failed.", ERR_ENTAP_INPUT_PARSE);
            }
        }

        // Test run DIAMOND if user input databases
        if (has_input(INPUT_FLAG_DATABASE)) {
            if (!ModDiamond::is_executable(dmnd_exe)) {
                throw ExceptionHandler("Databases have been selected for indexing. The test run of DIAMOND has failed!",
                                       ERR_ENTAP_INPUT_PARSE);
            }
        }
    }
    FS_dprint("Success!");
}

bool UserInput::has_input(ENTAP_INPUT_FLAGS input) {
     return !(mUserInputs[input].parsed_value).empty();
}


/**
 * ======================================================================
 * Function int get_supported_threads(boost::program_options::variables_map &user_map)
 *
 * Description          - Gets threads supported by system and compared with
 *                        user selection
 *                      - If thread support is lower, set to that
 *
 * Notes                - None
 *
 * @param user_map      - Boost parsed input
 *
 * @return              - thread number
 *
 * =====================================================================
 */
int UserInput::get_supported_threads() {

    uint32       supported_threads;
    int          threads;
    ent_input_uint_t          user_threads;

    supported_threads = std::thread::hardware_concurrency();
    user_threads = get_user_input<ent_input_uint_t >(INPUT_FLAG_THREADS);
    // assuming positive
    if ((uint32) user_threads > supported_threads) {
        FS_dprint("WARNING specified thread number is larger than available threads,"
                                        "setting threads to " + std::to_string(supported_threads));
        threads = supported_threads;
    } else {
        threads = user_threads;
    }
    return threads;
}

std::queue<char> UserInput::get_state_queue() {
    ent_input_str_t state_str;
    std::queue<char> out_queue;

    if (has_input(INPUT_FLAG_STATE)) {
        state_str = get_user_input<ent_input_str_t>(INPUT_FLAG_STATE);
        for (char c : state_str) {
            out_queue.push(c);
        }
    }
    return out_queue;   // check on return end if empty
}

std::string UserInput::get_target_species_str() {
    ent_input_str_t input_species;

    if (has_input(INPUT_FLAG_SPECIES)) {
        input_species = get_user_input<ent_input_str_t>(INPUT_FLAG_SPECIES);
        process_user_species(input_species);
        return input_species;
    } else return "";
}

vect_str_t UserInput::get_contaminants() {
    ent_input_multi_str_t output_contams;

    if (has_input(INPUT_FLAG_CONTAMINANT)) {
        output_contams = get_user_input<ent_input_multi_str_t>(INPUT_FLAG_CONTAMINANT);
        for (std::string &contam : output_contams) {
            if (contam.empty()) continue;
            process_user_species(contam);
        }
    }
    return output_contams;
}

vect_str_t UserInput::get_uninformative_vect() {
    ent_input_multi_str_t output_uninform;

    if (has_input(INPUT_FLAG_UNINFORMATIVE)) {
        output_uninform  = get_user_input<ent_input_multi_str_t>(INPUT_FLAG_UNINFORMATIVE);
    } else {
        return DEFAULT_UNINFORMATIVE;
    }

    for (std::string & str : output_uninform) {
        LOWERCASE(str);
    }
    return output_uninform;
}

std::string UserInput::get_user_transc_basename() {
    ent_input_str_t user_transcriptome;

    user_transcriptome = get_user_input<ent_input_str_t>(INPUT_FLAG_TRANSCRIPTOME);
    return mpFileSystem->get_filename(user_transcriptome, false);
}

std::vector<FileSystem::ENT_FILE_TYPES> UserInput::get_user_output_types() {
    std::vector<FileSystem::ENT_FILE_TYPES> ret;

    for (ent_input_uint_t val :get_user_input<ent_input_multi_int_t >(INPUT_FLAG_OUTPUT_FORMAT)) {
        ret.push_back(static_cast<FileSystem::ENT_FILE_TYPES>(val));
    }
    return ret;
}

std::string UserInput::getBIN_PATH_DEFAULT() {
    return BIN_PATH_DEFAULT;
}

std::string UserInput::getDATABASE_DIR_DEFAULT() {
    return DATABASE_DIR_DEFAULT;
}

const std::string &UserInput::getEGG_SQL_DB_FILENAME() {
    return EGG_SQL_DB_FILENAME;
}

const std::string &UserInput::getEGG_DMND_FILENAME() {
    return EGG_DMND_FILENAME;
}

const std::string &UserInput::getEGG_SQL_DB_DEFAULT() {
    return EGG_SQL_DB_DEFAULT;
}

const std::string &UserInput::getEGG_DMND_DEFAULT() {
    return EGG_DMND_DEFAULT;
}

const std::string &UserInput::getENTAP_DATABASE_BIN_DEFAULT() {
    return ENTAP_DATABASE_BIN_DEFAULT;
}

const std::string &UserInput::getENTAP_DATABASE_SQL_DEFAULT() {
    return ENTAP_DATABASE_SQL_DEFAULT;
}

// Returns false if unable to determine whether we want to run frame selection
bool UserInput::run_frame_selection(QueryData *queryData, bool &run_frame_selection) {
    FS_dprint("Determining if we want to run frame selection...");
    bool blastp;    // User input runP/blastp (TRUE) or runN (false) or config (false)
    bool ret;

    ret = true;
    blastp = has_input(INPUT_FLAG_RUNPROTEIN);

    if (queryData == nullptr) {
        FS_dprint("ERROR Unable to determine, nullptr!");
        ret = false;
    } else{

        if (blastp && queryData->is_protein_data()) {
            FS_dprint("NO Protein sequences input AND runP, skipping frame selection");
            run_frame_selection = false;
        } else if (!blastp) {
            FS_dprint("NO Blastx/runN selected, skipping frame selection");
            run_frame_selection = false;
        } else {
            FS_dprint("YES, run frame selection");
            run_frame_selection = true;
        }
    }

    return ret;
}

bool UserInput::run_expression_filtering() {
    FS_dprint("Determining if we want to run expression analysis");

    if (has_input(INPUT_FLAG_ALIGN)) {
        FS_dprint("YES, alignment file input from user");
        return true;
    } else {
        FS_dprint("NO, no alignment file specified from user");
        return false;
    }
}

ent_input_str_t UserInput::get_entap_database_path(EntapDatabase::DATABASE_TYPE type) {
    ent_input_str_t ret;

    switch (type) {
        case EntapDatabase::ENTAP_SERIALIZED:
            ret = get_user_input<ent_input_str_t>(INPUT_FLAG_ENTAP_DB_BIN);
            break;

        case EntapDatabase::ENTAP_SQL:
            ret = get_user_input<ent_input_str_t>(INPUT_FLAG_ENTAP_DB_SQL);
            break;

        default:
            // UNKNOWN database type
            break;
    }
    return ret;
}
