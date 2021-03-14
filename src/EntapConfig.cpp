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


//*********************** Includes *****************************
#include "EntapConfig.h"
#include "database/EggnogDatabase.h"
#include "TerminalCommands.h"
#include "database/BuscoDatabase.h"
//**************************************************************

namespace entapConfig {

    // ***************** Local Variables ***************************
    EntapDatabase           *pEntapDatabase;  // Pointer to EnTAP database (binary or sql)
    std::string              binDir;          // Absolute path to binary directory
    std::string              dataDir;         // Absolute path to database directory
    std::string              rootDir;         // Absolute path to EnTAP outfiles directory
    FileSystem              *pFileSystem;     // Pointer to EnTAP filesystem
    UserInput               *pUserInput;      // Pointer to User Data
    // *************************************************************

    // ***************** Local Prototype Functions******************
    void init_entap_database();
    void init_diamond_index(std::string &diamond_exe, uint16 threads, vect_str_t &compiled_databases);
    void init_eggnog(uint16 threads, std::string &diamond_exe);
    void init_busco();
    // *************************************************************


    /**
     * ======================================================================
     * Function void init_entap(UserInput *input, FileSystem *filesystem)
     *
     * Description          - Entry into configurating EnTAP
     *                      - Responsible for downloading EnTAP databases (taxonomic,
     *                        Gene Ontology, UniProt), DIAMOND configuring, and EggNOG download
     *
     * Notes                - Entry
     *
     * @param input         - Pointer to flags from user
     * @param filesystem    - Pointer to EnTAP filesystem
     *
     * @return              - None
     *
     * =====================================================================
     */
    void execute_main(UserInput *input, FileSystem *filesystem) {

        uint16                             threads;             // supported threads for execution
        ent_input_str_t                    diamond_exe;         // DIAMOND executable path
        ent_input_multi_str_t              compiled_databases;  // databases input from user

        FS_dprint("Entering EnTAP Config");

        pUserInput = input;
        pFileSystem = filesystem;

        // Get directory to output databases to. Same as "outfiles" path
        rootDir = filesystem->get_root_path();

        binDir  = PATHS(rootDir, input->getBIN_PATH_DEFAULT());
        dataDir = PATHS(rootDir, input->getDATABASE_DIR_DEFAULT());

        // Create paths
        pFileSystem->create_dir(rootDir);    // Should have already been created, confirm
        pFileSystem->create_dir(binDir);
        pFileSystem->create_dir(dataDir);

        threads = (uint16)pUserInput->get_supported_threads();
        diamond_exe = pUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_DIAMOND_EXE);

        if (pUserInput->has_input(INPUT_FLAG_DATABASE)) {
            compiled_databases = pUserInput->get_user_input<ent_input_multi_str_t>(INPUT_FLAG_DATABASE);
        }

        try {

            init_entap_database();

            init_diamond_index(diamond_exe, threads, compiled_databases);

            init_eggnog(threads, diamond_exe);

            // Disable BUSCO config for now
//            init_busco();

        } catch (ExceptionHandler &e){
            SAFE_DELETE(pEntapDatabase);
            throw ExceptionHandler(e.what(), e.getErr_code());
        }
        SAFE_DELETE(pEntapDatabase);
        FS_dprint("Configuration complete!");
    }

    /**
     * ======================================================================
     * Function init_diamond_index(std::string diamond_exe,std::string out_path,
     *                             uint16 threads)
     *
     * Description          - Responsible for indexing user specified FASTA formatted
     *                        database for DIAMOND usage
     *
     * Notes                - Utilizes DIAMIND executable
     *
     * @param diamond_exe   - Path to DIAMOND exe
     * @param out_path      - Database out directory
     * @param threads       - Thread number
     *
     * @return              - None
     *
     * =====================================================================
     */
    void init_diamond_index(std::string &diamond_exe, uint16 threads, vect_str_t &compiled_databases) {
        FS_dprint("Preparing to index database(s) with Diamond...");

        std::string indexed_path;       // Absolute path to final, DMND indexed database
        std::string std_out;            // Standard output (out, err) from execution
        std::string index_command;      // DMND indexing command
        std::stringstream log_msg;      // Message to print to EnTAP log file

        if (compiled_databases.empty()) {
            FS_dprint("No databases selected, skipping");
            // Exit routine if no databases need to be configured!!
            return;
        }

        pFileSystem->format_stat_stream(log_msg, "DIAMOND Database Configuration");

        for (std::string &fasta_path: compiled_databases) {
            TerminalData terminalData = TerminalData();

            indexed_path = PATHS(binDir,pFileSystem->get_filename(fasta_path, false));
            std_out      = indexed_path + FileSystem::EXT_STD;


            // TODO change for updated databases
            // Does database already exist?
            if (pFileSystem->file_exists(indexed_path + FileSystem::EXT_DMND)) {

                FS_dprint("File found at " + indexed_path + ".dmnd, skipping...");
                log_msg << "DIAMOND database skipped, exists at: " << indexed_path << std::endl;

            } else {
                // NO, database does not exist
                // Clear log
                pFileSystem->delete_file(std_out + FileSystem::EXT_ERR);
                pFileSystem->delete_file(std_out + FileSystem::EXT_OUT);

                index_command =
                        diamond_exe + " makedb --in " + fasta_path +
                        " -d "      + indexed_path +
                        " -p "      +std::to_string(threads);

                terminalData.command       = index_command;
                terminalData.base_std_path = std_out;
                terminalData.print_files   = true;
                terminalData.suppress_std_err = false;

                if (TC_execute_cmd(terminalData) != 0) {
                    throw ExceptionHandler("Error indexing database at: " + fasta_path + "\nDIAMOND Error: " + terminalData.err_stream,
                                           ERR_ENTAP_INIT_INDX_DATABASE);
                }
                FS_dprint("Database successfully indexed to: " + indexed_path + FileSystem::EXT_DMND);
                log_msg << "DIAMOND database generated to: " << indexed_path << FileSystem::EXT_DMND << std::endl;
            }
        } // END FOR

        std::string temp = log_msg.str();
        pFileSystem->print_stats(temp);
    }


    /**
     * ======================================================================
     * Function init_eggnog(uint16 threads, std::string &diamond_exe)
     *
     * Description          - Downloads and configured required EggNOG databases
     *                        (SQL database and FASTA database)
     *                      - Configured FASTA database for DIAMOND if necessary
     *
     * Notes                - Utilizes DIAMOND executable
     *
     * @param threads       - Thread number for execution
     * @param diamond_exe   - Method to execute DIAMOND
     *
     * @return              - None
     *
     * =====================================================================
     */
    void init_eggnog(uint16 threads, std::string &diamond_exe) {
        std::string sql_outpath;                // Absolute path to EggNOG sql output
        std::string fasta_outpath;              // Absolute path to EggNOG fasta output
        std::string dmnd_outpath;               // Absolute path to converted FASTA -> DMND
        std::string user_egg_dmnd;              // User input EggNOG DMND database
        std::string user_egg_sql;               // User input EggNOG SQL database
        std::string err_msg;                    // Error mMessage from execution
        std::string index_cmd;                  // DMND indexing command
        std::string std_out;                    // Standard output (err, out) from execution
        std::stringstream log_msg;              // Message to print to EnTAP log file

        FS_dprint("Ensuring EggNOG databases exist...");

        pFileSystem->format_stat_stream(log_msg, "EggNOG Database Configuration");

        // Generate database to allow downloading
        EggnogDatabase eggnogDatabase = EggnogDatabase(pFileSystem, pEntapDatabase, nullptr);

#if EGGNOG_MAPPER
        std::string eggnog_cmd;

        eggnog_cmd =
                "python " + eggnog_exe +
                " none -y";

        if (!pFileSystem->file_exists(eggnog_exe)) {
            throw ExceptionHandler("Eggnog download path does not exist at: " +
                                   eggnog_exe, ERR_ENTAP_INIT_EGGNOG);
        }
        if (TC_execute_cmd(eggnog_cmd) != 0) {
            throw ExceptionHandler("EggNOG command: " + eggnog_cmd,ERR_ENTAP_INIT_EGGNOG);
        }
        FS_dprint("Success! EggNOG databases verified");
        return;
#endif
        // setup outpath
        sql_outpath   = PATHS(dataDir, pUserInput->getEGG_SQL_DB_FILENAME());
        fasta_outpath = PATHS(pFileSystem->get_temp_outdir(), "eggnog_fasta_temp.fa");
        dmnd_outpath  = PATHS(binDir, pUserInput->getEGG_DMND_FILENAME());
        user_egg_dmnd = pUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_DMND_DB);
        user_egg_sql  = pUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_EGG_SQL_DB);

        // Check if SQL database already exists
        if (!pFileSystem->file_exists(user_egg_sql) && !pFileSystem->file_exists(sql_outpath)) {
            // No, path does not. Download it
            if (eggnogDatabase.download(EggnogDatabase::EGGNOG_SQL, sql_outpath) != EggnogDatabase::ERR_EGG_OK) {
                // Error in download
                err_msg = "Unable to download EggNOG sql database from FTP to: " +
                        sql_outpath + "\nError: " + eggnogDatabase.print_err();
                throw ExceptionHandler(err_msg,ERR_ENTAP_INIT_EGGNOG);
            } else {
                // Downloaded successfully
                FS_dprint("Success! EggNOG SQL database downloaded to: " + sql_outpath);
                log_msg << "EggNOG SQL database written to: " + sql_outpath << std::endl;
            }
        } else {
            // Already exists, skip
            std::string path;
            if (pFileSystem->file_exists(user_egg_sql)) path = user_egg_sql;
            if (pFileSystem->file_exists(sql_outpath)) path = sql_outpath;
            FS_dprint("EggNOG SQL database already exists at: " + path +
                " skipping");
            log_msg << "EggNOG SQL Database skipped, exists at: " << path << std::endl;
        }

        // Check if DIAMOND EggNOG database exists
        if (!pFileSystem->file_exists(user_egg_dmnd) && !pFileSystem->file_exists(dmnd_outpath)) {
            // No, does not exist, need to generate from FASTA
            if (eggnogDatabase.download(EggnogDatabase::EGGNOG_FASTA, fasta_outpath) != EggnogDatabase::ERR_EGG_OK) {
                // Error in download
                err_msg = "Unable to get EggNOG FASTA from FTP to: " +
                          fasta_outpath + "\nError: " + eggnogDatabase.print_err();
                throw ExceptionHandler(err_msg,ERR_ENTAP_INIT_EGGNOG);
            }

            // Now, index for DIAMOND
            FS_dprint("Success! EggNOG FASTA downloaded, indexing for DIAMOND...");
            log_msg << "EggNOG FASTA database written to: " + fasta_outpath << std::endl;

            // TODO move to DIAMOND module
            index_cmd =
                    diamond_exe + " makedb --in " + fasta_outpath +
                    " -d "      + dmnd_outpath +
                    " -p "      +std::to_string(threads);
            std_out = fasta_outpath + FileSystem::EXT_STD;

            TerminalData terminalData = TerminalData();

            terminalData.command       = index_cmd;
            terminalData.base_std_path = std_out;
            terminalData.print_files   = true;
            terminalData.suppress_std_err = false;

            if (TC_execute_cmd(terminalData) != 0) {
                throw ExceptionHandler("Error indexing database at: " + dmnd_outpath + "\nError:" +
                                       terminalData.err_stream, ERR_ENTAP_INIT_EGGNOG);
            }
            FS_dprint("Success! EggNOG DIAMOND database indexed at: " + dmnd_outpath);
            log_msg << "DIAMOND EggNOG database written to: " + dmnd_outpath << std::endl;

        } else {
            // Yes, file already exists
            std::string path;
            if (pFileSystem->file_exists(user_egg_dmnd)) path = user_egg_dmnd;
            if (pFileSystem->file_exists(dmnd_outpath)) path = dmnd_outpath;
            FS_dprint("EggNOG DIAMOND database already exists at: " + path);
            log_msg << "EggNOG DIAMOND database skipped, exists at: " << path << std::endl;
        }
        // Print to log/debug
        FS_dprint("Success! All EggNOG files verified");
        std::string temp = log_msg.str();
        pFileSystem->print_stats(temp);
    }

    /**
     * ======================================================================
     * Function init_entap_database()
     *
     * Description          - Downloads and configured required EggNOG databases
     *                        (SQL database and FASTA database)
     *                      - Configured FASTA database for DIAMOND if necessary
     *
     * Notes                - Requires Internet connection
     *
     *
     * @return              - None
     *
     * =====================================================================
     */
    void init_entap_database() {
        bool              generate_databases;    // Whether user would like to generate rather than download
        vect_uint16_t     databases;             // User defined types of databases to configure (SQL or Bin)
        std::string       config_outpath;        // Path to check against (from config file)
        std::string       database_outpath;      // Path to print to (also acts as default path)
        std::stringstream log_msg;               // Log mMessage to print to EnTAP log file
        EntapDatabase::DATABASE_TYPE database_type; // Type of database to configure/download (EntapDatabase.h)
        EntapDatabase::DATABASE_ERR database_err;   // Database error (EntapDatabase.h)

        FS_dprint("Initializing EnTAP database...");
        pFileSystem->format_stat_stream(log_msg, "EnTAP Database Configuration");

        pEntapDatabase = new EntapDatabase(pFileSystem);
        if (pEntapDatabase == nullptr) {
            throw ExceptionHandler("Unable to allocate Entap Database memory", ERR_ENTAP_MEM_ALLOC);
        }

        // If user would like to generate databases rather than download them from ftp(default)
        generate_databases = pUserInput->has_input(INPUT_FLAG_DATABASE_GENERATE);

        // Check which databases they want (will always have this input, default = 0)
        databases = pUserInput->get_user_input<ent_input_multi_int_t>(INPUT_FLAG_DATABASE_TYPE);

        // Download or generate databases
        FS_dprint("Beginning to download/generate databases...");
        for (uint16 data : databases) {
            // Check database
            database_type = static_cast<EntapDatabase::DATABASE_TYPE>(data);

            // Set outpaths and paths to check against
            switch (database_type) {
                case EntapDatabase::ENTAP_SERIALIZED:
                    FS_dprint("Generating/downloading Serialized database...");
                    config_outpath  = pUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_ENTAP_DB_BIN);
                    database_outpath = PATHS(rootDir, pUserInput->getENTAP_DATABASE_BIN_DEFAULT());
                    break;

                case EntapDatabase::ENTAP_SQL:
                    FS_dprint("Generating/downloading SQL database");
                    config_outpath   = pUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_ENTAP_DB_SQL);;
                    database_outpath = PATHS(rootDir, pUserInput->getENTAP_DATABASE_SQL_DEFAULT());
                    break;

                default:
                    FS_dprint("WARNING Unrecognized database code: " + std::to_string(data));
                    continue;
            }

            // First check if this database already exists in the config outpath or default outpath
            if (pFileSystem->file_exists(config_outpath) || pFileSystem->file_exists(database_outpath)) {
                std::string path;
                if (pFileSystem->file_exists(config_outpath)) path = config_outpath;
                if (pFileSystem->file_exists(database_outpath)) path = database_outpath;
                FS_dprint("File already exists at: " + path);
                log_msg << "Database skipped, already exists at: " << path << std::endl;
                continue; // Don't redownload
            }

            // Need to generate/download file!
            if (generate_databases) {
                FS_dprint("EntapConfig: Generating database to: " + database_outpath + "...");
                database_err = pEntapDatabase->generate_database(database_type, database_outpath);
            } else {
                FS_dprint("EntapConfig: Downloading database to: " + database_outpath);
                database_err = pEntapDatabase->download_database(database_type, database_outpath);
            }

            // Check if successful
            if (database_err == EntapDatabase::ERR_DATA_OK) {
                FS_dprint("Success! Database written to: " + database_outpath);
                log_msg << "Database written to: " + database_outpath << std::endl;
            } else {
                // Fatal if any databases fail
                throw ExceptionHandler(pEntapDatabase->print_error_log(), ERR_ENTAP_INIT_DATA_GENERIC);
            }
        } // END FOR

        // print to log
        std::string temp = log_msg.str();
        pFileSystem->print_stats(temp);
        SAFE_DELETE(pEntapDatabase);
    }

    /**
     * ======================================================================
     * Function init_busco()
     *
     * Description          - Downloads and configures required files for the
     *                        BUSCO module from user input (OrthoDB datasets)
     *
     * Notes                - Require Internet connection
     *
     *
     * @return              - None
     *
     * =====================================================================
     */
    void init_busco() {
        std::string busco_outpath;          // Absolute path to download BUSCO database to
        ent_input_multi_str_t busco_databases;         // BUSCO database tag/URL input by user
        std::string busco_out_filename;     // Filename of output BUSCO database
        BuscoDatabase::BUSCO_DB_ERR busco_db_err;   // BUSCO database error code
        std::stringstream log_msg;               // Log message to print to EnTAP log file

        if (!pUserInput->has_input(INPUT_FLAG_BUSCO_DATABASE)) {
            FS_dprint("No BUSCO database input by the user, skipping...");
            return; // EXIT - user has not input a BUSCO database
        }

        FS_dprint("Initializing BUSCO database...");
        busco_db_err = BuscoDatabase::ERR_DATA_OK;
        pFileSystem->format_stat_stream(log_msg, "BUSCO Database Configuration");
        busco_databases = pUserInput->get_user_input<ent_input_multi_str_t>(INPUT_FLAG_BUSCO_DATABASE);

        for (std::string &busco_database : busco_databases) {
            if (pFileSystem->is_url(busco_database)) {
                busco_out_filename = pFileSystem->get_filename(busco_database, false);
            } else {
                busco_out_filename = busco_database;
            }
            busco_outpath = PATHS(dataDir, busco_out_filename);

            if (pFileSystem->file_exists(busco_outpath)) {
                FS_dprint("BUSCO database already exists at: " + busco_outpath);
                log_msg << "BUSCO Database skipped, already exists at: " << busco_outpath << std::endl;

            } else {
                // Database does not exist, download it
                BuscoDatabase buscoDatabase = BuscoDatabase(pFileSystem);
                busco_db_err = buscoDatabase.download_database(busco_database, busco_outpath);

                if (busco_db_err == BuscoDatabase::ERR_DATA_OK) {
                    FS_dprint("Success! BUSCO database downloaded to: " + busco_outpath);
                    log_msg << "BUSCO Database written to: " << busco_outpath << std::endl;
                } else {
                    // Fatal if database download fails
                    throw ExceptionHandler(buscoDatabase.print_error_log(), ERR_ENTAP_INIT_BUSCO_GENERIC);
                }
            }
        }

        std::string temp = log_msg.str();
        pFileSystem->print_stats(temp);
    }
}
