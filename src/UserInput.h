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

#ifndef ENTAP_USERINPUT_H
#define ENTAP_USERINPUT_H

//*********************** Includes *****************************
#include <boost/program_options/variables_map.hpp>
#include "EntapGlobals.h"

//**************************************************************

//*********************** Defines ******************************
#define DESC_COMING_SOON    "Coming soon!"
#define DESC_HELP           "Print all the help options for this version of EnTAP!"
#define DESC_CONFIG         "Configure EnTAP for execution later.\n"                    \
                            "If this is your first time running EnTAP run this first!"  \
                            "This will perform the following:\n"                        \
                            "    - Downloading EnTAP taxonomic database\n"              \
                            "    - Downloading Gene Ontology term database\n"           \
                            "    - Formatting any database you would like for diamond"
#define DESC_RUN_PROTEIN    "Execute EnTAP functionality through blastp\n"              \
                            "Note, if your input sequences are nucleotide, they will be"\
                            "frame selected automatically."
#define DESC_RUN_NUCLEO     "Execute EnTAP functionality through blastx\n"              \
                            "This will not frame select your sequences and will run them"\
                            "through each stage of the pipeline as nucelotide sequences"
#define DESC_INTER_DATA     "Select which databases you would like for InterProScan"    \
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
#define DESC_ONTOLOGY_FLAG  "Specify the ontology software you would like to use\n"     \
                            "Note: it is possible to specify more than one! Just use"   \
                            "multiple --ontology flags\n"                               \
                            "Specify flags as follows:\n"                               \
                            "    1. EggNOG (default)\n"                                 \
                            "    2. InterProScan"
#define DESC_GRAPHING       "Check whether or not your system supports graphing.\n"     \
                            "This option does not require any other flags and will"     \
                            "just check whether the version of Python being used has"   \
                            "MatPlotLib accessible."
#define DESC_OUT_FLAG       "Specify the output directory you would like the data to"   \
                            " be saved to."
#define DESC_DATABASE       "Provide the paths to the databases you would like to use\n"\
                            "For running: ensure the databases selected are .dmnd"      \
                            "formatted.\n"                                              \
                            "For configuration: ensure the databases are FASTA format\n"\
                            "Note: if your databases are not NCBI or Uniprot\n"         \
                            "databases, taxonomic filtering might not be able to pull"  \
                            "the species information!"
#define DESC_ONT_LEVELS     "Specify the Gene Ontology levels you would like printed\n" \
                            "Default: 0, 3, 4\n"                                        \
                            "A level of 0 means that every term will be printed!"       \
                            "It is possible to specify multiple flags as well with\n"   \
                            "multiple --level flags\n"                                  \
                            "Example: --level 0 --level 3 --level 1"
#define DESC_FPKM           "Specify the FPKM threshold with expression analysis\n"     \
                            "EnTAP will filter out transcripts below this value!"
#define DESC_EVAL           "Specify the E-Value that will be used as a cutoff during"  \
                            "similarity searching"
#define DESC_THREADS        "Specify the number of threads that will be used throughout\n"
#define DESC_SINGLE_END     "Specify this flag if your BAM/SAM file was generated\n"    \
                            "through single-end reads\n"                                \
                            "Note: this is only required in expression analysis\n"      \
                            "Default: paired-end"
#define DESC_ALIGN_FILE     "Specify the path to the BAM/SAM file for expression analysis"
#define DESC_CONTAMINANT    "Specify the contaminants you would like to filter out"     \
                            "from similarity searching\n"                               \
                            "Note: since hits are based upon a multitide of factors"    \
                            "a contaminant might be the best hit for a query!\n"        \
                            "Contaminants can be selected by species (homo_sapiens)"    \
                            "or through a specific taxon (homo)\n"                      \
                            "If your taxon is more than one word just replace the"      \
                            "spaces with underscores (_)"
#define DESC_TRIM           "Trim the input sequences to the first space\n"             \
                            "This may help with readability later on with TSV files\n"  \
                            "Example:\n"                                                \
                            ">TRINITY_231.1 Protein Information\n"                      \
                            "will become...\n"                                          \
                            ">TRINITY_231.1\n"
#define DESC_QCOVERAGE      "Select the minimum query coverage to be allowed during"    \
                            "similarity searching"
#define DESC_TCOVERAGE      "Select the minimum target coverage to be allowed during"   \
                            "similarity searching"
#define DESC_EXE_PATHS      "Specify path to the entap_config.txt file that will"       \
                            "be used to find all of the executables!"
#define DESC_DATA_OUT       "Specify the outpath of databases formatted during"         \
                            "configuration\n"                                           \
                            "Note: only used in configuration stage"
#define DESC_TAXON          "Specify the type of species/taxon you are analyzing and"   \
                            "would like hits closer in taxonomic relevance to be"       \
                            "favored (based on NCBI Taxonomic Database)\n"              \
                            "Note: formatting works just like with the contaminants"
#define DESC_STATE          "Specify the state of execution (EXPERIMENTAL)\n"           \
                            "More information is available in the documentation\n"      \
                            "This flag may have undesired affects and may not run properly!"
#define DESC_INPUT_TRAN     "Path to the input transcriptome file"
#define DESC_COMPLET_PROT   "Select this option if all of your sequences are complete"  \
                            "proteins.\n"                                               \
                            "At this point, this option will merely flag the sequences"
#define DESC_OVERWRITE      "Select this option if you would like to overwrite previous"\
                            "files\n"                                                   \
                            "Note: do NOT use this if you would like to pickup from"    \
                            "a previous run!"
#define DESC_UNINFORMATIVE  "Path to a list of keywords that should be used to specify" \
                            " uninformativeness of hits during similarity searching. "  \
                            "Generally something along the lines of 'hypothetical' or " \
                            "'unknown' are used. Each term should be on a new line of " \
                            "the file being linked to.\nExample (defaults):\n"          \
                            "    -conserved\n"                                          \
                            "    -predicted\n"                                          \
                            "    -unknown\n"                                            \
                            "    -hypothetical\n"                                       \
                            "    -putative\n"                                           \
                            "    -unidentified\n"                                       \
                            "    -uncultured\n"                                         \
                            "    -uninformative\n"                                      \
                            "Without the extra spaces and hyphen. EnTAP will take the " \
                            "each line as a new uninformative word!"
#define DESC_NOCHECK        "Use this flag if you don't want your input to EnTAP verifed."\
                            " This is not advised to use! Your run may fail later on "  \
                            "if inputs are not checked"
//**************************************************************

//*********************** Typedefs/Enum ************************
typedef std::vector<std::string> databases_t;

enum SPECIES_FLAGS {
    SPECIES,
    CONTAMINANT
};
//******************** Prototype Functions *********************
boost::program_options::variables_map parse_arguments_boost(int, const char**);
bool verify_user_input(boost::program_options::variables_map&);
void print_user_input(boost::program_options::variables_map &map, std::string&, std::string&);
bool check_key(std::string&);
std::unordered_map<std::string,std::string> parse_config(std::string&,std::string&);
void generate_config(std::string&);
void verify_databases(boost::program_options::variables_map&);
void verify_species (boost::program_options::variables_map&, SPECIES_FLAGS);
void init_exe_paths(std::unordered_map<std::string, std::string> &, std::string);
std::string get_exe_path(boostPO::variables_map&);
bool verify_interpro(std::string);
void process_user_species(std::string&);
void verify_uninformative(std::string&);
void verify_state(std::string&, bool, std::vector<uint16>&);
std::pair<bool,std::string> verify_software(uint8&, std::vector<uint16>&);

//**************************************************************



//*********************** Constants ****************************
const std::string INTER_TIGR            = "tigrfam";
const std::string INTER_SFLD            = "sfld";
const std::string INTER_PRODOM          = "prodom";
const std::string INTER_HAMAP           = "hamap";
const std::string INTER_PFAM            = "pfam";
const std::string INTER_SMART           = "smart";
const std::string INTER_CDD             = "cdd";
const std::string INTER_PROSITE_PROF    = "prositeprofiles";
const std::string INTER_PROSITE_PAT     = "prositepatterns";
const std::string INTER_SUPERFAMILY     = "superfamily";
const std::string INTER_PRINTS          = "prints";
const std::string INTER_PANTHER         = "panther";
const std::string INTER_GENE            = "gene3d";
const std::string INTER_PIRSF           = "pirsf";
const std::string INTER_COILS           = "coils";
const std::string INTER_MOBI            = "mobidblite";

const fp32 DEFAULT_QCOVERAGE               = 50.0;
const fp32 DEFAULT_TCOVERAGE               = 50.0;
const fp32 COVERAGE_MIN                    = 0.0;
const fp32 COVERAGE_MAX                    = 100.0;
const fp32 E_VALUE                         = 1e-5;
const fp32 RSEM_FPKM_DEFAULT               = 0.5;
const fp32 FPKM_MIN                        = 0.0;
const fp32 FPKM_MAX                        = 100.0;
const uint8 MAX_DATABASE_SIZE              = 5;
const std::string DEFAULT_STATE            = "+";
const std::string INTERPRO_DEFAULT         = INTER_PFAM;
const std::string OUTFILE_DEFAULT          = PATHS(boostFS::current_path(),"outfiles");
const std::string BAM_EXT                  = ".bam";
const std::string SAM_EXT                  = ".sam";


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
// All paths based around main EnTAP directory
namespace Defaults {
    const std::string RSEM_DEFAULT_EXE         = "/libs/RSEM-1.3.0/";   // Directory
    const std::string GENEMARK_DEFAULT_EXE     = "/libs/gmst_linux_64/gmst.pl";
    const std::string DIAMOND_DEFAULT_EXE      = "/libs/diamond-0.8.31/bin/diamond";
    const std::string EGG_EMAPPER_DEFAULT      = "/libs/eggnog-mapper/emapper.py";
    const std::string EGG_DOWNLOAD_DEFAULT     = "/libs/eggnog-mapper/download_eggnog_data.py";
    const std::string EGG_SQL_DB_DEFAULT       = "/libs/eggnog-mapper/data/eggnog.db";
    const std::string INTERPRO_DEF_EXE         = "/libs/interproscan-5.22-61.0/interproscan.sh";
    const std::string TAX_DOWNLOAD_DEF         = "/src/download_tax.py";
    const std::string GRAPH_SCRIPT_DEF         = "/src/entap_graphing.py";
}

//**************************************************************






#endif //ENTAP_USERINPUT_H
