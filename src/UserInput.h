/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
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

#ifdef USE_BOOST
#include <boost/program_options/options_description.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options/variables_map.hpp>
#else
#include <tclap/CmdLine.h>
#include <boost/any.hpp>          // Include any boost library with tclap
#endif

//**************************************************************


namespace Defaults {
    const std::string RSEM_DEFAULT_EXE         = "/libs/RSEM-1.3.0/";   // Directory
    const std::string GENEMARK_DEFAULT_EXE     = "/libs/gmst_linux_64/gmst.pl";
    const std::string DIAMOND_DEFAULT_EXE      = "/libs/diamond-0.8.31/bin/diamond";
    const std::string EGG_SQL_DB_FILENAME      = "eggnog.db";
    const std::string EGG_DMND_FILENAME        = "eggnog_proteins.dmnd";
    const std::string INTERPRO_DEF_EXE         = "interproscan.sh";
    const std::string TAX_DOWNLOAD_DEF         = "/src/download_tax.py";
    const std::string GRAPH_SCRIPT_DEF         = "/src/entap_graphing.py";
    const std::string BIN_PATH_DEFAULT         = "/bin";
    const std::string DATABASE_DIR_DEFAULT     = "/databases";
    const std::string ENTAP_DATABASE_SQL_FILENAME    = "entap_database.db";
    const std::string ENTAP_DATABASE_SQL_GZ    = "entap_database.db.gz";
    const std::string ENTAP_DATABASE_SERIAL_FILENAME = "entap_database.bin";
    const std::string ENTAP_DATABASE_SERIAL_GZ = "entap_database.bin.gz";
    const std::string ENTAP_DATABASE_BIN_DEFAULT = PATHS(BIN_PATH_DEFAULT, ENTAP_DATABASE_SERIAL_FILENAME);
    const std::string ENTAP_DATABASE_SQL_DEFAULT = PATHS(DATABASE_DIR_DEFAULT, ENTAP_DATABASE_SQL_FILENAME);
    const std::string EGG_SQL_DB_DEFAULT         = PATHS(DATABASE_DIR_DEFAULT, EGG_SQL_DB_FILENAME);
    const std::string EGG_DMND_DEFAULT           = PATHS(BIN_PATH_DEFAULT, EGG_DMND_FILENAME);
}

class UserInput {

public:
    UserInput(int argc, const char** argv);
    ~UserInput();

    bool has_input(const std::string&);
    pair_str_t get_config_path();
    void set_pFileSystem(FileSystem *_pFileSystem);
    std::unordered_map<std::string,std::string> parse_config(pair_str_t&);
    bool verify_user_input();
    int get_supported_threads();
    std::queue<char> get_state_queue();
    std::string get_target_species_str();
    vect_str_t get_contaminants();
    vect_str_t get_uninformative_vect();
    std::string get_user_transc_basename();

    template<class T>
    T get_user_input(const std::string &key) {
#ifdef USE_BOOST
        if (_user_inputs.count(key)) {
            return _user_inputs[key].as<T>();
        } else {
            return T();
        }
#else // Use TCLAP
        if (has_input(key)) {
            return boost::any_cast<T>(_user_inputs[key]);
        } else {
            return T();
        }
#endif
    }

    const std::string INPUT_FLAG_TAG           = "out-dir";
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
    const std::string INPUT_FLAG_SPECIES       = "taxon";
    const std::string INPUT_FLAG_QCOVERAGE     = "qcoverage";
    const std::string INPUT_FLAG_TCOVERAGE     = "tcoverage";
    const std::string INPUT_FLAG_COMPLETE      = "complete";
    const std::string INPUT_FLAG_GO_LEVELS     = "level";
    const std::string INPUT_FLAG_EXE_PATH      = "paths";
    const std::string INPUT_FLAG_FPKM          = "fpkm";
    const std::string INPUT_FLAG_CONTAM        = "contam";
    const std::string INPUT_FLAG_E_VAL         = "e";
    const std::string INPUT_FLAG_HELP          = "help";
    const std::string INPUT_FLAG_VERSION       = "version";
    const std::string INPUT_FLAG_TRANSCRIPTOME = "input";
    const std::string INPUT_FLAG_DATABASE      = "database";
    const std::string INPUT_FLAG_GRAPH         = "graph";
    const std::string INPUT_FLAG_TRIM          = "trim";
    const std::string INPUT_FLAG_STATE         = "state";
    const std::string INPUT_FLAG_SINGLE_END    = "single-end";
    const std::string INPUT_FLAG_THREADS       = "threads";
    const std::string INPUT_FLAG_UNINFORM      = "uninformative";
    const std::string INPUT_FLAG_NOCHECK       = "no-check";
    const std::string INPUT_FLAG_GENERATE      = "data-generate";
    const std::string INPUT_FLAG_DATABASE_TYPE = "data-type";

private:
    enum SPECIES_FLAGS {
        SPECIES,
        CONTAMINANT
    };

#ifdef USE_BOOST
    void parse_arguments_boost(int, const char**);
#else
    void parse_arguments_tclap(int, const char **);
#endif
    void print_user_input();
    bool check_key(std::string&);
    void generate_config(std::string&);
    void verify_databases(bool);
    void verify_species (SPECIES_FLAGS, EntapDatabase*);
    void init_exe_paths(std::unordered_map<std::string, std::string> &, std::string);
    void process_user_species(std::string&);
    void verify_uninformative(std::string&);
    void verify_state(std::string&, bool, std::vector<uint16>&);
    std::pair<bool,std::string> verify_software(uint8&, std::vector<uint16>&);
    std::string get_executable_dir();

#ifdef USE_BOOST
    boostPO::variables_map _user_inputs;
#else
    std::map<std::string, boost::any> _user_inputs;     // Header only library for any map
#endif

    const fp32 DEFAULT_QCOVERAGE               = 50.0;
    const fp32 DEFAULT_TCOVERAGE               = 50.0;
    const fp32 COVERAGE_MIN                    = 0.0;
    const fp32 COVERAGE_MAX                    = 100.0;
    const fp64 E_VALUE                         = 1e-5;
    const uint32 DEFAULT_THREADS               = 1;
    const fp32 RSEM_FPKM_DEFAULT               = 0.5;
    const fp32 FPKM_MIN                        = 0.0;
    const fp32 FPKM_MAX                        = 100.0;
    const uint8 MAX_DATABASE_SIZE              = 5;
    const std::string DEFAULT_STATE            = "+";
    const std::string OUTFILE_DEFAULT          = PATHS(FileSystem::get_cur_dir(),"outfiles");

    // Enter as lowercase
    const std::vector<std::string> INFORMATIVENESS {
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
    const std::string KEY_EGGNOG_SQL_DB        = "eggnog_sql_database";
    const std::string KEY_EGGNOG_DMND          = "eggnog_dmnd_database";
    const std::string KEY_ENTAP_DATABASE_BIN   = "entap_database_bin_path";
    const std::string KEY_ENTAP_DATABASE_SQL   = "entap_database_sql_path";
    const std::string KEY_GRAPH_SCRIPT         = "entap_graphing_script";

    FileSystem *_pFileSystem;
    bool        _is_config;

//**************************************************************
};





#endif //ENTAP_USERINPUT_H
