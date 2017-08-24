//
// Created by harta55 on 2/1/17.
//

#include <map>
#include "EntapConfig.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include "pstream.h"
#include "boost/filesystem.hpp"
#include "EntapGlobals.h"
#include "ExceptionHandler.h"
#include "EntapExecute.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/program_options/variables_map.hpp>
#include <thread>
#include <csv.h>
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"

namespace boostAR = boost::archive;

namespace entapConfig {

    enum InitStates {
        INIT_TAX            = 0x01,
        INIT_UNIPROT        = 0x02,
        INIT_NCBI           = 0x04,
        INIT_DATABASE       = 0x08,
        INIT_DIAMOND_INDX   = 0x16,
        INIT_EXIT           = 0x32
    };

    InitStates               state;
    std::vector<std::string> _compiled_databases;
    std::string              _bin_dir;
    std::string              _data_dir;
    std::string              _outpath;
    std::string              _cur_dir;

    void init_entap(boost::program_options::variables_map user_map, std::string exe_path) {

        std::string                        database_outdir;
        int                                threads; // change
        std::pair<std::string,std::string> exe_pair;
        std::vector<std::string>           ncbi_vect;
        std::vector<std::string>           uniprot_vect;
        std::vector<std::string>           database_vect;


        if (user_map.count(ENTAP_CONFIG::INPUT_FLAG_DATA_OUT)) {
            database_outdir = user_map[ENTAP_CONFIG::INPUT_FLAG_DATA_OUT].as<std::string>();
        } else database_outdir = exe_path;

        _cur_dir  = boostFS::path(boost::filesystem::current_path()).string();
        _outpath  = (boostFS::path(_cur_dir) / user_map["tag"].as<std::string>()).string();
        _bin_dir  = PATHS(database_outdir, ENTAP_CONFIG::BIN_PATH);
        _data_dir = PATHS(database_outdir, ENTAP_CONFIG::DATABASE_DIR);

        if (database_outdir.empty()) database_outdir = PATHS(_cur_dir,ENTAP_CONFIG_DIR);

        boostFS::create_directories(database_outdir);
        boostFS::create_directories(_bin_dir);
        boostFS::create_directories(_data_dir);

        threads = get_supported_threads(user_map);

        ncbi_vect = user_map[ENTAP_CONFIG::INPUT_FLAG_NCBI_1].as<std::vector<std::string>>();
        uniprot_vect = user_map[ENTAP_CONFIG::INPUT_FLAG_UNIPROT].as<std::vector<std::string>>();
        if (user_map.count("database")) {
            database_vect = user_map["database"].as<std::vector<std::string>>();
            _compiled_databases = database_vect;
        }

        // while state != EXIT_STATE
        try {
            init_taxonomic(database_outdir);
            init_go_db(database_outdir,_data_dir);
#if NCBI_UNIPROT
            init_uniprot(uniprot_vect, exe_path);
            init_ncbi(ncbi_vect,exe_path);
#endif
            init_diamond_index(DIAMOND_EXE,database_outdir, threads);
            init_eggnog(EGG_DOWNLOAD_EXE);

        }catch (ExceptionHandler &e) {
            throw ExceptionHandler(e.what(), e.getErr_code());
        }
    }

    void init_taxonomic(std::string &exe) {
        print_debug("Downloading taxonomic database...");
        //TODO Integrate gzip/zlib

        std::string tax_bin;
        std::string tax_path;
        std::string tax_command;
        std::unordered_map<std::string, std::string> tax_data_map;

        tax_path  = (boostFS::path(exe) / boostFS::path(TAX_DATABASE_PATH)).string();

        if (file_exists(TAX_DB_PATH)) {
            tax_bin = TAX_DB_PATH;
            print_debug("Tax database binary found at: " + tax_bin + " skipping...");
            return;
            // TODO update!
        } else {
            print_debug("Tax database not found at: " + tax_bin + " downloading...");
            tax_command = "perl " + TAX_DOWNLOAD_EXE + " " + tax_path;
            if (execute_cmd(tax_command) != 0) {
                throw ExceptionHandler("Command: " + tax_command, ENTAP_ERR::E_INIT_TAX_DOWN);
            }
            print_debug("Success! File written to " + tax_path);
        }

        print_debug("Indexing taxonomic database...");
        std::ifstream infile(tax_path);
        std::string line;
        try {
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                std::string lineage, sci_name, tax_id;
                std::getline(iss, sci_name, '\t');
                std::getline(iss, tax_id, '\t');
                std::getline(iss, lineage, '\t');

                std::string id_lineage = tax_id+"||"+lineage;
                tax_data_map.emplace(sci_name,id_lineage);
            }
        } catch (std::exception &e) {
            throw ExceptionHandler(e.what(), ENTAP_ERR::E_INIT_TAX_INDEX);
        }

        print_debug("Success! Writing file to "+ tax_bin);
        try{
            {
                std::ofstream ofs(tax_bin);
                boostAR::binary_oarchive oa(ofs);
                oa << tax_data_map;
            }
        } catch (std::exception &e) {
            throw ExceptionHandler(e.what(), ENTAP_ERR::E_INIT_TAX_SERIAL);
        }
        print_debug("Success!");
    }

    void init_go_db(std::string &exe, std::string database_path) {
        print_debug("Initializing GO terms database...");

        std::string go_db_path;
        std::string go_term_path;
        std::string go_graph_path;
        std::string go_database_zip;
        std::string go_database_out;

        if (file_exists(GO_DB_PATH)) {
            print_debug("Database found at: " + GO_DB_PATH + " skipping creation");
            return;
        } else {
            print_debug("Database NOT found at: " + GO_DB_PATH + " , downloading...");
            go_db_path = PATHS(exe, ENTAP_CONFIG::GO_DB_PATH_DEF);
        }
        go_database_zip = (boostFS::path(database_path) / boostFS::path(GO_DATA_NAME)).string();
        go_database_out = (boostFS::path(database_path) / boostFS::path(GO_DIR)).string();
        if (file_exists(go_database_zip)) {
            boostFS::remove(go_database_zip);
        }
        try {
            download_file(GO_DATABASE_FTP,go_database_zip);
            decompress_file(go_database_zip,database_path,1);
            boostFS::remove(go_database_zip);
        } catch (ExceptionHandler const &e ){ throw e;}

        go_term_path = (boostFS::path(go_database_out) / boostFS::path(GO_TERM_FILE)).string();
        go_graph_path = (boostFS::path(go_database_out) / boostFS::path(GO_GRAPH_FILE)).string();
        if (!file_exists(go_term_path) || !file_exists(go_graph_path)) {
            throw ExceptionHandler("GO term files must be at: " + go_term_path + " and " + go_graph_path +
                                   " in order to configure database", ENTAP_ERR::E_INIT_GO_SETUP);
        }

        io::CSVReader<6, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(go_graph_path);
        std::string index,root,branch, temp, distance, temp2;
        std::map<std::string,std::string> distance_map;
        const std::string BIOLOGICAL = "6679", MOLECULAR = "2892", CELLULAR = "311";
        while (in.read_row(index,root,branch, temp, distance, temp2)) {
            if (root.compare(BIOLOGICAL) == 0 || root.compare(MOLECULAR) == 0 || root.compare(CELLULAR) ==0) {
                if (distance_map.find(branch) == distance_map.end()) {
                    distance_map.emplace(branch,distance);
                } else {
                    if (distance.empty()) continue;
                    long cur = std::stoi(distance_map[branch]);
                    long query = std::stoi(distance);
                    if (query > cur) distance_map[branch] = distance;
                }
            }
        }
        std::map<std::string,struct_go_term> go_map;
        std::string num,term,cat,go,ex,ex1,ex2;
        io::CSVReader<7, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in2(go_term_path);
        while (in2.read_row(num,term,cat,go,ex,ex1,ex2)) {
            std::string lvl = distance_map[num];
            struct_go_term go_data;
            go_data.category = cat; go_data.level=lvl;go_data.term=term;
            go_data.go_id=go;
            go_map[go] = go_data;
        }
        boostFS::remove_all(go_database_out);
        print_debug("Success! Writing file to "+ go_db_path + "...");
        try{
            {
                std::ofstream ofs(go_db_path);
                boostAR::binary_oarchive oa(ofs);
                oa << go_map;
            }
        } catch (std::exception &e) {
            throw ExceptionHandler(e.what(), ENTAP_ERR::E_INIT_GO_SETUP);
        }
        print_debug("Success!");
    }

    // may handle differently than ncbi with formatting
    void init_uniprot(std::vector<std::string> &flags, std::string exe) {
        // TODO setup go term/interpro... integration, date tag, use bool input
        print_debug("Parsing uniprot databases...");
        if (flags.empty()) {
            print_debug("No Uniprot databases selected");
            return;
        }
        std::string ftp_address;
        std::string uniprot_bin = exe + "/" + ENTAP_CONFIG::BIN_PATH + "uniprot_";
        std::string uniprot_data = exe + ENTAP_CONFIG::UNIPROT_BASE_PATH;

        for (auto &flag : flags) {
            if (flag == ENTAP_CONFIG::INPUT_UNIPROT_NULL) return;
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
        std::string ncbi_data = exe + ENTAP_CONFIG::NCBI_BASE_PATH;
        for (auto &flag : flags) {
            if (flag == ENTAP_CONFIG::INPUT_UNIPROT_NULL) return;
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

    void init_diamond_index(std::string diamond_exe,std::string out_path,int threads) {
        print_debug("Preparing to index database(s) with Diamond...");

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
            boostFS::remove(std_out+".err");
            boostFS::remove(std_out+".out");

            // TODO change for updated databases
            if (file_exists(indexed_path + ".dmnd")) {
                print_debug("File found at " + indexed_path + ".dmnd, skipping...");
                continue;
            }
            index_command =
                    diamond_exe + " makedb --in " + item +
                    " -d "      + indexed_path +
                    " -p "      +std::to_string(threads);

            print_debug("Executing DIAMOND command:\n" + index_command);
            if (execute_cmd(index_command,std_out) != 0) {
                throw ExceptionHandler("Error indexing database at: " + item,
                                       ENTAP_ERR::E_INIT_INDX_DATABASE);
            }
            print_debug("Database successfully indexed to: " + indexed_path + ".dmnd");
        }
    }

    void init_eggnog(std::string eggnog_exe) {

        std::string eggnog_cmd;

        eggnog_cmd =
                "python " + eggnog_exe +
                " none -y";

        if (!file_exists(eggnog_exe)) {
            throw ExceptionHandler("Eggnog download path does not exist at: " +
                                   eggnog_exe, ENTAP_ERR::E_INIT_EGGNOG);
        }
        print_debug("Executing eggnog download...\n" + eggnog_cmd);
        if (execute_cmd(eggnog_cmd) != 0) {
            throw ExceptionHandler("Error in executing eggnog download",ENTAP_ERR::E_INIT_EGGNOG);
        }
    }

    std::string download_file(std::string flag, std::string &path, std::string temp) {
        // TODO libcurl

        std::string ftp_address;
        int status;
        std::string download_command;
        std::string output_path;

        output_path = path + flag + ".gz";

        if (flag == ENTAP_CONFIG::INPUT_UNIPROT_SWISS) {
            ftp_address = UNIPROT_FTP_SWISS;

        } else {
            throw ExceptionHandler("Invalid uniprot flag", ENTAP_ERR::E_INPUT_PARSE);
        }

        download_command = "wget -O "+ output_path + " " + ftp_address;
        print_debug("Downloading uniprot: " + flag + " database from " +
                  ftp_address + "...");
        status = execute_cmd(download_command);
        if (status != 0) {
            throw ExceptionHandler("Error in downloading uniprot database", ENTAP_ERR::E_INIT_TAX_DOWN);
        }
        print_debug("File successfully downloaded to: " + output_path);
        return output_path;
    }

    std::string download_file(const std::string &ftp, std::string &out_path) {
        boostFS::path path(out_path);
        std::string download_command = "wget -O "+ out_path + " " + ftp;
        print_debug("Downloading through wget: file from " + ftp + "...");
        execute_cmd(download_command);
        print_debug("Success! File printed to: " + out_path);
        return out_path;
    }

    void decompress_file(std::string file_path, std::string out_path,short flag) {
        print_debug("Decompressing file at: " + file_path);

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
                    file_path, ENTAP_ERR::E_INIT_TAX_DOWN);
        }
        print_debug("File at: " + file_path + " successfully decompressed");
    }

    int update_database(std::string file_path) {
        return 0;
    }
}
