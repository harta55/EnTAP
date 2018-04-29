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


//*********************** Includes *****************************
#include "EntapConfig.h"
//**************************************************************

namespace entapConfig {

    enum InitStates {
        INIT                = 0,
        INIT_ENTAP_DATABASE,
        INIT_DIAMOND_INDX  ,
        INIT_EGGNOG        ,
        INIT_EXIT          ,
    };

    InitStates               state;
    std::vector<std::string> _compiled_databases;
    EntapDatabase            *_pEntapDatabase;
    std::string              _bin_dir;
    std::string              _data_dir;
    std::string              _outpath;
    FileSystem               *_pFileSystem;
    UserInput                *_pUserInput;

    //-----------------------FTP PATHS---------------------------//
    const std::string UNIPROT_FTP_SWISS = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";
    const std::string UNIPROT_FTP_TREMBL = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz";

    //****************** Local Prototype Functions******************
    void init_entap_database(std::string&);
    void init_uniprot(std::vector<std::string>&, std::string);
    void init_ncbi(std::vector<std::string>&, std::string);
    void init_diamond_index(std::string,std::string,int);
    int update_database(std::string);
    void init_eggnog(std::string);
    void handle_state(void);


    /**
     * ======================================================================
     * Function void init_entap(boost::program_options::variables_map user_map,
     *                          std::string exe_path)
     *
     * Description          - Entry into configurating EnTAP
     *                      - Responsible for downloading EnTAP databases (taxonomic,
     *                        Gene Ontology), DIAMOND configuring, and EggNOG download
     *
     * Notes                - Entry
     *
     * @param user_map      - Boost parsed user input flags
     * @param exe_path      - Path to EnTAP executable and main directory
     *
     * @return              - None
     *
     * =====================================================================
     */
    void execute_main(UserInput *input, FileSystem *filesystem) {

        std::string                        database_outdir;
        int                                threads; // change
        std::pair<std::string,std::string> exe_pair;
        std::vector<std::string>           ncbi_vect;
        std::vector<std::string>           uniprot_vect;
        std::vector<std::string>           database_vect;

        _pUserInput = input;
        _pFileSystem = filesystem;

        // Get directory to output databases to. Same as "outfiles" path
        database_outdir = filesystem->get_root_path();

        _bin_dir  = PATHS(database_outdir, Defaults::BIN_PATH_DEFAULT);
        _data_dir = PATHS(database_outdir, Defaults::DATABASE_DIR_DEFAULT);

        _pFileSystem->create_dir(database_outdir);
        _pFileSystem->create_dir(_bin_dir);
        _pFileSystem->create_dir(_data_dir);

        threads = _pUserInput->get_supported_threads();

        if (_pUserInput->has_input(UInput::INPUT_FLAG_DATABASE)) {
            database_vect = _pUserInput->get_user_input<vect_str_t>(UInput::INPUT_FLAG_DATABASE);
            _compiled_databases = database_vect;
        }

        // while state != EXIT_STATE
        while (state != INIT_EXIT) {
            try {
                switch (state) {
                    case INIT_ENTAP_DATABASE:
                        init_entap_database(database_outdir);
                        break;
#if NCBI_UNIPROT
                    case INIT_UNIPROT:
                        init_uniprot(uniprot_vect, exe_path);
                        break;
                    case INIT_NCBI:
                        iknit_ncbi(ncbi_vect,exe_path);
                        break;
#endif
                    case INIT_DIAMOND_INDX:
                        init_diamond_index(DIAMOND_EXE, database_outdir, threads);
                        break;
                    case INIT_EGGNOG:
                        init_eggnog(EGG_DOWNLOAD_EXE);
                        break;
                    default:
                        break;
                }
            handle_state();
            } catch (ExceptionHandler &e) {
                throw ExceptionHandler(e.what(), e.getErr_code());
            }
        }

        SAFE_DELETE(_pEntapDatabase);
        FS_dprint("Configuration complete!");
    }


#if NCBI_UNIPROT
    // may handle differently than ncbi with formatting
    void init_uniprot(std::vector<std::string> &flags, std::string exe) {
        // TODO setup go term/interpro... integration, date tag, use bool input
        print_debug("Parsing uniprot databases...");
        if (flags.empty()) {
            print_debug("No Uniprot databases selected");
            return;
        }
        std::string ftp_address;
        std::string uniprot_bin = exe + "/" + UInput::BIN_PATH + "uniprot_";
        std::string uniprot_data = exe + UInput::UNIPROT_BASE_PATH;

        for (auto &flag : flags) {
            if (flag == UInput::INPUT_UNIPROT_NULL) return;
            std::string diamond_path = uniprot_bin + flag + ".dmnd";
            std::string database_path = uniprot_data + flag + ".fasta";
            if (file_exists(database_path)) {
                print_debug("Database at: " + database_path + " found, updating...");
                update_database(database_path);
                _compiled_databases.push_back(database_path);
            } else {
                print_debug("Database at: " + database_path + " not found, downloading...");
                try {
                    std::string temp_path = download_file(flag, database_path);
                    decompress_file(temp_path,temp_path,0);
                    _compiled_databases.push_back(database_path);
                } catch (ExceptionHandler &e) {throw e;}
            }
        }
    }

    void init_ncbi(std::vector<std::string> &flags, std::string exe) {
        // TODO setup go term/interpro... integration, date tag, use bool input
        print_debug("Parsing NCBI databases...");
        if (flags.empty()) {
            print_debug("No NCBI databases selected");
            return;
        }
        std::string ftp_address;
        std::string ncbi_data = exe + UInput::NCBI_BASE_PATH;
        for (auto &flag : flags) {
            if (flag == UInput::INPUT_UNIPROT_NULL) return;
            std::string database_path = ncbi_data + flag + ".fasta";
            if (file_exists(database_path)) {
                print_debug("Database at: " + database_path + " found, updating...");
                update_database(database_path);
                _compiled_databases.push_back(database_path);
            } else {
                print_debug("Database at: " + database_path + " not found, downloading...");
                try {
                    std::string temp_path = download_file(flag, database_path);
                    decompress_file(temp_path,temp_path,0);
                    _compiled_databases.push_back(database_path);
                } catch (ExceptionHandler &e) {throw e;}
            }
        }
    }
#endif


    /**
     * ======================================================================
     * Function init_diamond_index(std::string diamond_exe,std::string out_path,
     *                             int threads)
     *
     * Description          - Responsible for indexing user specified FASTA formatted
     *                        database for DIAMOND usage
     *
     * Notes                - Utilizes script in /src
     *
     * @param diamond_exe   - Path to DIAMOND exe
     * @param out_path      - Database out directory
     * @param threads       - Thread number
     *
     * @return              - None
     *
     * =====================================================================
     */
    void init_diamond_index(std::string diamond_exe,std::string out_path,int threads) {
        FS_dprint("Preparing to index database(s) with Diamond...");

        std::string filename;
        std::string indexed_path;
        std::string std_out;
        std::string index_command;

        if (_compiled_databases.empty()) return;

        for (std::string item : _compiled_databases) {
            boostFS::path path(item);
            filename     = path.filename().stem().string();
            indexed_path = PATHS(out_path,filename);
            std_out      = indexed_path + "_index";
            _pFileSystem->delete_file(std_out + FileSystem::EXT_ERR);
            _pFileSystem->delete_file(std_out + FileSystem::EXT_OUT);

            // TODO change for updated databases
            if (_pFileSystem->file_exists(indexed_path + ".dmnd")) {
                FS_dprint("File found at " + indexed_path + ".dmnd, skipping...");
                continue;
            }
            index_command =
                    diamond_exe + " makedb --in " + item +
                    " -d "      + indexed_path +
                    " -p "      +std::to_string(threads);

            FS_dprint("Executing DIAMOND command:\n" + index_command);
            if (TC_execute_cmd(index_command, std_out) != 0) {
                throw ExceptionHandler("Error indexing database at: " + item,
                                       ERR_ENTAP_INIT_INDX_DATABASE);
            }
            FS_dprint("Database successfully indexed to: " + indexed_path + ".dmnd");
        }
    }


    /**
     * ======================================================================
     * Function init_eggnog(std::string eggnog_exe)
     *
     * Description          - Ensure EggNOG databases are downloaded
     *
     * Notes                - Only will download DIAMOND database
     *
     * @param eggnog_exe    - Path to EggNOG download python script
     *
     * @return              - None
     *
     * =====================================================================
     */
    void init_eggnog(std::string eggnog_exe) {
        FS_dprint("Ensuring EggNOG-mapper databases exist...");

        std::string eggnog_cmd;

        eggnog_cmd =
                "python " + eggnog_exe +
                " none -y";

        if (!_pFileSystem->file_exists(eggnog_exe)) {
            throw ExceptionHandler("Eggnog download path does not exist at: " +
                                   eggnog_exe, ERR_ENTAP_INIT_EGGNOG);
        }
        if (TC_execute_cmd(eggnog_cmd) != 0) {
            throw ExceptionHandler("EggNOG command: " + eggnog_cmd,ERR_ENTAP_INIT_EGGNOG);
        }
        FS_dprint("Success! EggNOG databases verified");
    }

    // TODO update databases
    int update_database(std::string file_path) {
        return 0;
    }

    void handle_state(void) {
        state = static_cast<InitStates>(state+1);
    }

    void init_entap_database(std::string &out_dir) {
        bool generate_databases;    // Whether user would like to generate rather tahn download
        vect_uint16_t databases;
        std::string config_outpath;    // Path to check against (from config file)
        std::string database_outpath;  // Path to print to (also acts as default path)
        EntapDatabase::DATABASE_TYPE database_type;
        EntapDatabase::DATABASE_ERR database_err;

        FS_dprint("Initializing EnTAP database...");

        _pEntapDatabase = new EntapDatabase(_pFileSystem);
        if (_pEntapDatabase == nullptr) {
            throw ExceptionHandler("Unable to allocate Entap Database memory", ERR_ENTAP_MEM_ALLOC);
        }

        // If user would like to generate databases rather than download them from ftp(default)
        generate_databases = _pUserInput->has_input(UInput::INPUT_FLAG_GENERATE);

        // Check which databases they want (will always have this input, default = 0)
        databases = _pUserInput->get_user_input<vect_uint16_t>(UInput::INPUT_FLAG_DATABASE_TYPE);

        // Download or generate databases
        FS_dprint("Beginning to download/generate databases...");
        for (uint16 data : databases) {
            // Check database
            database_type = static_cast<EntapDatabase::DATABASE_TYPE>(data);

            // Set outpaths and paths to check against
            switch (database_type) {
                case EntapDatabase::ENTAP_SERIALIZED:
                    FS_dprint("Generating/downloading Serialized database...");
                    config_outpath  = ENTAP_DATABASE_BIN_PATH;
                    database_outpath = PATHS(out_dir, Defaults::ENTAP_DATABASE_BIN_DEFAULT);
                    break;

                case EntapDatabase::ENTAP_SQL:
                    FS_dprint("Generating/downloading SQL database");
                    config_outpath   = ENTAP_DATABASE_SQL_PATH;
                    database_outpath = PATHS(out_dir, Defaults::ENTAP_DATABASE_SQL_DEFAULT);
                    break;

                default:
                    FS_dprint("Unrecognized database code: " + std::to_string(data));
                    continue;
            }

            // First check if this database already exists in the config outpath or default outpath
            if (_pFileSystem->file_exists(config_outpath) || _pFileSystem->file_exists(database_outpath)) {
                FS_dprint("File already exists at: " + config_outpath);
                continue; // Don't redownload
            }

            // Need to generate/download file!
            if (generate_databases) {
                FS_dprint("EntapConfig: Generating database to: " + database_outpath + "...");
                database_err = _pEntapDatabase->generate_database(database_type, database_outpath);
            } else {
                FS_dprint("EntapConfig: Downloading database to: " + database_outpath);
                database_err = _pEntapDatabase->download_database(database_type, database_outpath);
            }

            // Check if successful
            if (database_err == EntapDatabase::ERR_DATA_OK) {
                FS_dprint("Success! Database written to: " + database_outpath);
            } else {
                throw ExceptionHandler("Error in getting database " + std::to_string(data) +
                    ". Database Error: " + _pEntapDatabase->print_error_log(database_err), ERR_ENTAP_INIT_DATA_GENERIC);
            }
        }
    }
}
