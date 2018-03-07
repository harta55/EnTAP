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


//*********************** Includes *****************************
#include "EntapConfig.h"
#include <sys/stat.h>
#include <unordered_map>
#include "pstream.h"
#include "boost/filesystem.hpp"
#include "EntapGlobals.h"
#include "ExceptionHandler.h"
#include "EntapExecute.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/program_options/variables_map.hpp>
#include <thread>
#include <csv.h>
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "FileSystem.h"
//**************************************************************

namespace entapConfig {

    enum InitStates {
        INIT                = 0,
        INIT_TAX           ,
        INIT_UNIPROT       ,
        INIT_NCBI          ,
        INIT_DATABASE      ,
        INIT_DIAMOND_INDX  ,
        INIT_EGGNOG        ,
        INIT_GO            ,
        INIT_EXIT          ,
    };

    InitStates               state;
    std::vector<std::string> _compiled_databases;
    std::string              _bin_dir;
    std::string              _data_dir;
    std::string              _outpath;
    FileSystem               *_pFileSystem;
    UserInput                *_pUserInput;

    //-----------------------FTP PATHS---------------------------//
    const std::string GO_DATABASE_FTP =
            "http://archive.geneontology.org/latest-full/go_monthly-termdb-tables.tar.gz";
    const std::string UNIPROT_FTP_SWISS = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";
    const std::string UNIPROT_FTP_TREMBL = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz";
    const std::string GO_TERM_FILE             = "term.txt";
    const std::string GO_GRAPH_FILE            = "graph_path.txt";
    const std::string GO_DATA_NAME             = "go_monthly-termdb-tables.tar.gz";
    const std::string GO_DIR                   = "go_monthly-termdb-tables/";

    // ID's used for GO level determination
    const std::string BIOLOGICAL_LVL = "6679";
    const std::string MOLECULAR_LVL  = "2892";
    const std::string CELLULAR_LVL   = "311";


    //****************** Local Prototype Functions******************
    void init_taxonomic(std::string&);
    void init_uniprot(std::vector<std::string>&, std::string);
    void init_ncbi(std::vector<std::string>&, std::string);
    void init_diamond_index(std::string,std::string,int);
    std::string download_file(std::string, std::string&,std::string&);
    std::string download_file(const std::string &,std::string&);
    std::string get_output_dir();
    void decompress_file(std::string,std::string,short);
    int update_database(std::string);
    void init_go_db(std::string&,std::string);
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

        // Get directory to output databases to. Databases may already exist though!
        // If no database_out flag set, then defaults to outfiles directory
        database_outdir = get_output_dir();

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
                    case INIT_TAX:
                        init_taxonomic(database_outdir);
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
                    case INIT_GO:
                        init_go_db(database_outdir, _data_dir);
                        break;
                    default:
                        break;
                }
            handle_state();
            } catch (ExceptionHandler &e) {
                throw ExceptionHandler(e.what(), e.getErr_code());
            }
        }
    }


    /**
     * ======================================================================
     * Function void init_taxonomic(std::string &exe)
     *
     * Description          - Responsible for downloading NCBI Taxonomic database
     *
     * Notes                - Utilizes Perl script in /src
     *
     * @param user_map      - Boost parsed user input flags
     *
     * @return              - None
     *
     * =====================================================================
     */
    void init_taxonomic(std::string &exe) {
        FS_dprint("Downloading taxonomic database...");
        //TODO Integrate gzip/zlib

        std::string tax_bin_path;
        std::string tax_txt_path_default;
        std::string tax_bin_path_default;
        std::string tax_command;
        tax_serial_map_t tax_data_map;

        // Tax txt default is not in the config file, so have to check default
        tax_txt_path_default  = PATHS(exe, Defaults::TAX_DATABASE_TXT_DEFAULT);
        // Pull tax bin path from config file paths
        tax_bin_path          = TAX_DB_PATH;
        // If it doesn't exist in config file path, will check default path
        tax_bin_path_default  = PATHS(exe, Defaults::TAX_DATABASE_BIN_DEFAULT);

        if (_pFileSystem->file_exists(TAX_DB_PATH)) {
            FS_dprint("Tax database binary found at: " + tax_bin_path + " skipping...");
            return;
            // TODO update database and check date of download
        } else {
            // No bin database located at config file path!!
            FS_dprint("Tax database binary not found at: " + TAX_DB_PATH + " checking non-binary...");
            // Check if bin database exists at the output directory
            if (_pFileSystem->file_exists(tax_bin_path_default)) {
                FS_dprint("Tax database binary found at: " + tax_bin_path_default + " skipping...");
                return;
            }
            // Bin database not found! Check if text version does not exist. Then download
            if (!_pFileSystem->file_exists(tax_txt_path_default)) {
                FS_dprint("Tax database not found at: " + tax_txt_path_default + " downloading...");
                tax_command = "python " + TAX_DOWNLOAD_EXE + " -o " + tax_txt_path_default;
                if (execute_cmd(tax_command) != 0) {
                    throw ExceptionHandler("Command: " + tax_command, ERR_ENTAP_INIT_TAX_DOWN);
                }
                FS_dprint("Success! File written to " + tax_txt_path_default);
            } else FS_dprint("Non-binary database found at: " + tax_txt_path_default + " indexing...");
        }

        // Text version of tax database should exist at this point (either existing or downloaded)
        FS_dprint("Indexing taxonomic database...");
        std::ifstream infile(tax_txt_path_default);
        std::string line;
        try {
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                std::string lineage, sci_name, tax_id;
                std::getline(iss, sci_name, '\t');
                std::getline(iss, tax_id, '\t');
                std::getline(iss, lineage, '\t');
                LOWERCASE(lineage);        // Both should already be lowercase
                LOWERCASE(sci_name);

                TaxEntry taxEntry = {
                        tax_id,
                        lineage
                };
                tax_data_map.emplace(sci_name,taxEntry);
            }
        } catch (std::exception &e) {
            throw ExceptionHandler(e.what(), ERR_ENTAP_INIT_TAX_INDEX);
        }

        FS_dprint("Success! Writing file to "+ tax_bin_path_default);
        try{
            {
                std::ofstream ofs(tax_bin_path);
                boostAR::binary_oarchive oa(ofs);
                oa << tax_data_map;
            }
        } catch (std::exception &e) {
            throw ExceptionHandler(e.what(), ERR_ENTAP_INIT_TAX_SERIAL);
        }
        // Update extern
        TAX_DB_PATH = tax_bin_path;
        FS_dprint("Success!");
    }


    /**
     * ======================================================================
     * Function void init_go_db(std::string &exe, std::string database_path)
     *
     * Description          - Responsible for downloading and indexing a
     *                        mapping of the Gene Ontology database
     *
     * Notes                - Utilizes script in /src
     *
     * @param user_map      - Boost parsed user input flags
     *
     * @return              - None
     *
     * =====================================================================
     */
    void init_go_db(std::string &exe, std::string database_path) {
        FS_dprint("Initializing GO terms database...");

        std::string go_db_path;
        std::string go_term_path;
        std::string go_graph_path;
        std::string go_database_zip;
        std::string go_database_out;
        std::string lvl;

        if (_pFileSystem->file_exists(GO_DB_PATH)) {
            FS_dprint("Database found at: " + GO_DB_PATH + " skipping creation");
            return;
        } else {
            FS_dprint("Database NOT found at: " + GO_DB_PATH + " , downloading...");
            go_db_path = PATHS(exe, Defaults::GO_DATABASE_BIN_DEFAULT);
        }
        go_database_zip = PATHS(database_path, GO_DATA_NAME);
        go_database_out = PATHS(database_path, GO_DIR);
        if (_pFileSystem->file_exists(go_database_zip)) boostFS::remove(go_database_zip);
        try {
            download_file(GO_DATABASE_FTP,go_database_zip);
            decompress_file(go_database_zip,database_path,1);
            boostFS::remove(go_database_zip);
        } catch (ExceptionHandler const &e ){ throw e;}

        go_term_path  = PATHS(go_database_out, GO_TERM_FILE);
        go_graph_path = PATHS(go_database_out, GO_GRAPH_FILE);
        if (!_pFileSystem->file_exists(go_term_path) || !_pFileSystem->file_exists(go_graph_path)) {
            throw ExceptionHandler("GO term files must be at: " + go_term_path + " and " + go_graph_path +
                                   " in order to configure database", ERR_ENTAP_INIT_GO_INDEX);
        }

        io::CSVReader<6, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(go_graph_path);
        std::string index,root,branch, temp, distance, temp2;
        std::map<std::string,std::string> distance_map;
        while (in.read_row(index,root,branch, temp, distance, temp2)) {
            if (root.compare(entapConfig::BIOLOGICAL_LVL) == 0     ||
                    root.compare(entapConfig::MOLECULAR_LVL) == 0  ||
                    root.compare(entapConfig::CELLULAR_LVL) ==0) {
                if (distance_map.find(branch) == distance_map.end()) {
                    distance_map.emplace(branch,distance);
                } else {
                    if (distance.empty()) continue;
                    fp32 cur   = std::stoi(distance_map[branch]);
                    fp32 query = std::stoi(distance);
                    if (query > cur) distance_map[branch] = distance;
                }
            }
        }
        go_serial_map_t go_map;
        std::string num,term,cat,go,ex,ex1,ex2;
        io::CSVReader<7, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in2(go_term_path);
        while (in2.read_row(num,term,cat,go,ex,ex1,ex2)) {
            lvl = distance_map[num];
            GoEntry go_data;
            go_data.category = cat; go_data.level=lvl;go_data.term=term;
            go_data.go_id=go;
            go_map[go] = go_data;
        }
        boostFS::remove_all(go_database_out);
        FS_dprint("Success! Writing file to "+ go_db_path + "...");
        try{
            {
                std::ofstream ofs(go_db_path);
                boostAR::binary_oarchive oa(ofs);
                oa << go_map;
            }
        } catch (std::exception &e) {
            throw ExceptionHandler(e.what(), ERR_ENTAP_INIT_GO_INDEX);
        }
        FS_dprint("Success!");
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
            boostFS::remove(std_out + FileSystem::EXT_ERR);
            boostFS::remove(std_out + FileSystem::EXT_OUT);

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
            if (execute_cmd(index_command,std_out) != 0) {
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

        std::string eggnog_cmd;

        eggnog_cmd =
                "python " + eggnog_exe +
                " none -y";

        if (!_pFileSystem->file_exists(eggnog_exe)) {
            throw ExceptionHandler("Eggnog download path does not exist at: " +
                                   eggnog_exe, ERR_ENTAP_INIT_EGGNOG);
        }
        FS_dprint("Executing eggnog download...\n" + eggnog_cmd);
        if (execute_cmd(eggnog_cmd) != 0) {
            throw ExceptionHandler("EggNOG command: " + eggnog_cmd,ERR_ENTAP_INIT_EGGNOG);
        }
    }

    std::string download_file(std::string flag, std::string &path, std::string temp) {
        // TODO libcurl

        std::string ftp_address;
        int status;
        std::string download_command;
        std::string output_path;

        output_path = path + flag + ".gz";

        download_command = "wget -O "+ output_path + " " + ftp_address;
        FS_dprint("Downloading uniprot: " + flag + " database from " +
                  ftp_address + "...");
        status = execute_cmd(download_command);
        if (status != 0) {
            throw ExceptionHandler("Error in downloading uniprot database", ERR_ENTAP_INIT_TAX_DOWN);
        }
        FS_dprint("File successfully downloaded to: " + output_path);
        return output_path;
    }

    std::string download_file(const std::string &ftp, std::string &out_path) {
        int status;

        boostFS::path path(out_path);
        std::string download_command = "wget -O "+ out_path + " " + ftp;
        FS_dprint("Downloading through wget: file from " + ftp + "...");
        status = execute_cmd(download_command);
        if (status != 0) {
            throw ExceptionHandler("Error in downloading " + ftp, ERR_ENTAP_INIT_GO_DOWNLOAD);
        }
        FS_dprint("Success! File printed to: " + out_path);
        return out_path;
    }


    /**
     * ======================================================================
     * Function decompress_file(std::string file_path, std::string out_path,short flag)
     *
     * Description          - Decompresses file using unix commands
     *
     * Notes                - Will be replaced with library
     *
     * @param file_path     - Path to file to be unzipped
     * @param out_path      - Output file path
     * @param flag          - Gzip/tar flag
     *
     * @return              - None
     *
     * =====================================================================
     */
    void decompress_file(std::string file_path, std::string out_path,short flag) {
        FS_dprint("Decompressing file at: " + file_path);

        std::string unzip_command;
        std::string std_out;
        int status;

        if (flag == 0){
            unzip_command = "gzip -d " + file_path;
        } else {
            unzip_command = "tar -xzf " + file_path + " -C " + out_path;
        }
        std_out = out_path + "_std";
        status = execute_cmd(unzip_command,std_out);
        if (status != 0) {
            throw ExceptionHandler("Error in unzipping database at " +
                    file_path, ERR_ENTAP_INIT_GO_UNZIP);
        }
        FS_dprint("File at: " + file_path + " successfully decompressed");
    }

    int update_database(std::string file_path) {
        return 0;
    }

    void handle_state(void) {
        state = static_cast<InitStates>(state+1);
    }

    std::string get_output_dir() {
        if (_pUserInput->has_input(UInput::INPUT_FLAG_DATA_OUT)) {
            return _pUserInput->get_user_input<std::string>(UInput::INPUT_FLAG_DATA_OUT);
        } return _pFileSystem->get_root_path();
    }
}
