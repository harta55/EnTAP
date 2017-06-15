//
// Created by harta55 on 2/1/17.
//

#include <map>
#include "EntapInit.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include "pstream.h"
#include "boost/filesystem.hpp"
#include "EntapConsts.h"
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

namespace boostFS = boost::filesystem;
namespace boostAR = boost::archive;
namespace Chrono = std::chrono;

namespace entapInit {

    enum InitStates {
        INIT_TAX            = 0x01,
        INIT_UNIPROT        = 0x02,
        INIT_NCBI           = 0x04,
        INIT_DATABASE       = 0x08,
        INIT_DIAMOND_INDX   = 0x16,
        INIT_EXIT           = 0x32
    };

    InitStates state;
    const boostFS::path current_path(boost::filesystem::current_path());
    std::list<std::string> _compiled_databases;

    void init_entap(boost::program_options::variables_map user_map, std::string exe_path,
            std::unordered_map<std::string,std::string> &config_map) {

        std::string outpath = current_path.string() + user_map["tag"].as<std::string>();
        boostFS::path bin_dir(exe_path + ENTAP_CONFIG::BIN_PATH);
        boostFS::path data_dir(exe_path + "/databases");
        bool bin_dir_state = (boostFS::create_directories(bin_dir));
        bool data_dir_state = (boostFS::create_directories(data_dir));

        int threads = get_supported_threads(user_map);

        std::vector<std::string> ncbi_vect, uniprot_vect, database_vect;
        ncbi_vect = user_map[ENTAP_CONFIG::INPUT_FLAG_NCBI_1].as<std::vector<std::string>>();
        uniprot_vect = user_map[ENTAP_CONFIG::INPUT_FLAG_UNIPROT].as<std::vector<std::string>>();
        if (user_map.count("database"))
            database_vect= user_map["database"].as<std::vector<std::string>>();

        std::string diamond_exe = entapExecute::init_exe_paths(config_map,exe_path);
        // while state != EXIT_STATE
        try {
            state = INIT_TAX;
            init_taxonomic(exe_path);
            init_go_db(exe_path,data_dir.string());
            init_uniprot(uniprot_vect, exe_path);
            init_ncbi(ncbi_vect,exe_path);
            init_diamond_index(diamond_exe,exe_path, threads);

        }catch (ExceptionHandler &e) {
            throw ExceptionHandler(e.what(), e.getErr_code());
        }
    }

    bool file_exists(const std::string &name) {
        struct stat buff;
        return (stat(name.c_str(), &buff) == 0);
    }

    void init_taxonomic(std::string &exe) {
        print_msg("Downloading taxonomic database...");
        //TODO Integrate gzip/zlib
        std::string tax_path = exe + ENTAP_CONFIG::TAX_DATABASE_PATH;

        if (!file_exists(tax_path)) {
            std::string tax_command = "perl " + exe + ENTAP_CONFIG::TAX_SCRIPT_PATH;
            redi::ipstream in(tax_command);
            in.close();
            int status = in.rdbuf()->status();
            if (status != 0) {
                std::cerr << "Error in downloading taxonomic database" << std::endl;
                throw ExceptionHandler("Error in downloading taxonomic database", ENTAP_ERR::E_INIT_TAX_DOWN);
            }
            print_msg("Success! File written to " + tax_path);
        } else {
            print_msg("Database found. Updating...");
            // TODO Update taxonomic database
            return;
        }

        print_msg("Indexing taxonomic database...");

        std::unordered_map<std::string, std::string> tax_data_map;
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

        print_msg("Success!");
        std::string tax_bin = exe + ENTAP_CONFIG::TAX_BIN_PATH;
        print_msg("Writing file to "+ tax_bin);
        // write map
        try{
            {
                std::ofstream ofs(tax_bin);
                boostAR::binary_oarchive oa(ofs);
                oa << tax_data_map;
            }
        } catch (std::exception &e) {
            throw ExceptionHandler(e.what(), ENTAP_ERR::E_INIT_TAX_SERIAL);
        }
        print_msg("Success!");

    }

    void init_go_db(std::string &exe, std::string database_path) {
        print_msg("Initializing GO terms database...");
        std::string go_db_path = exe + ENTAP_CONFIG::GO_DB_PATH;
        if (file_exists(go_db_path)) {
            print_msg("Database found at: " + go_db_path + " skipping creation");
            return;
        }
        const std::string GO_TERM_FILE = "term.txt";
        const std::string GO_GRAPH_FILE = "graph_path.txt";
        const std::string GO_DATABASE_FTP =
                "http://archive.geneontology.org/latest-full/go_monthly-termdb-tables.tar.gz";
        const std::string GO_DATA_NAME = "go_monthly-termdb-tables.tar.gz";
        const std::string GO_DIR = "go_monthly-termdb-tables/";
        std::string go_database_zip = database_path + "/" + GO_DATA_NAME;
        std::string go_database_out = database_path + "/" + GO_DIR;
        if (file_exists(go_database_zip)) {
            boostFS::remove(go_database_zip);
        }
        try {
            download_file(GO_DATABASE_FTP,go_database_zip);
            decompress_file(go_database_zip,database_path,1);
            boostFS::remove(go_database_zip);
        } catch (ExceptionHandler const &e ){ throw e;}

        std::string go_term_path = go_database_out + GO_TERM_FILE;
        std::string go_graph_path = go_database_out + GO_GRAPH_FILE;
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
        print_msg("Success!");
        print_msg("Writing file to "+ go_db_path);
        try{
            {
                std::ofstream ofs(go_db_path);
                boostAR::binary_oarchive oa(ofs);
                oa << go_map;
            }
        } catch (std::exception &e) {
            throw ExceptionHandler(e.what(), ENTAP_ERR::E_INIT_GO_SETUP);
        }
        print_msg("Success!");
    }

    // may handle differently than ncbi with formatting
    void init_uniprot(std::vector<std::string> &flags, std::string exe) {
        // TODO setup go term/interpro... integration, date tag, use bool input
        print_msg("Parsing uniprot databases...");
        if (flags.empty()) {
            print_msg("No Uniprot databases selected");
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
                print_msg("Database at: " + database_path + " found, updating...");
                update_database(database_path);
                _compiled_databases.push_back(database_path);
            } else {
                print_msg("Database at: " + database_path + " not found, downloading...");
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
        print_msg("Parsing NCBI databases...");
        if (flags.empty()) {
            print_msg("No NCBI databases selected");
            return;
        }
        std::string ftp_address;
        std::string ncbi_data = exe + ENTAP_CONFIG::NCBI_BASE_PATH;
        for (auto &flag : flags) {
            if (flag == ENTAP_CONFIG::INPUT_UNIPROT_NULL) return;
            std::string database_path = ncbi_data + flag + ".fasta";
            if (file_exists(database_path)) {
                print_msg("Database at: " + database_path + " found, updating...");
                update_database(database_path);
                _compiled_databases.push_back(database_path);
            } else {
                print_msg("Database at: " + database_path + " not found, downloading...");
                try {
                    std::string temp_path = download_file(flag, database_path);
                    decompress_file(temp_path,temp_path,0);
                    _compiled_databases.push_back(database_path);
                } catch (ExceptionHandler &e) {throw e;}
            }
        }
    }

    void init_diamond_index(std::string diamond_exe,std::string exe_path,int threads) {
        print_msg("Preparing to index database(s) with Diamond...");
        if (_compiled_databases.empty()) return;
        std::string bin_path = exe_path + "/" + ENTAP_CONFIG::BIN_PATH;
        for (std::string item : _compiled_databases) {
            boostFS::path path(item);
            std::string filename = path.filename().stem().string();
            std::string indexed_path = bin_path + filename;
            std::string std_out = bin_path + filename + "_index";
            boostFS::remove(std_out+".err");
            boostFS::remove(std_out+".out");

            // TODO change for updated databases
            if (file_exists(indexed_path + ".dmnd")) {
                print_msg("File found at " + indexed_path + ".dmnd, skipping...");
                continue;
            }
            std::string index_command = diamond_exe + " makedb --in " +
                item + " -d " + indexed_path + " -p "+std::to_string(threads);
            if (execute_cmd(index_command,std_out) != 0) {
                throw ExceptionHandler("Error indexing database at: " + item,
                                       ENTAP_ERR::E_INIT_INDX_DATABASE);
            }
            print_msg("Database successfully indexed to: " + indexed_path + ".dmnd");
        }
    }

    void print_msg(std::string msg) {
        Chrono::time_point<Chrono::system_clock> current = Chrono::system_clock::now();
        std::time_t time = Chrono::system_clock::to_time_t(current);
        std::string out_time(std::ctime(&time));
        std::ofstream log_file("debug.txt", std::ios::out | std::ios::app);
        log_file << out_time.substr(0,out_time.length()-1) << ": " + msg << std::endl;
        log_file.close();
    }

    int execute_cmd(std::string cmd, std::string out_path) {
        std::ofstream out_file(out_path+".out", std::ios::out | std::ios::app);
        std::ofstream err_file(out_path+".err", std::ios::out | std::ios::app);
        const redi::pstreams::pmode mode = redi::pstreams::pstdout|redi::pstreams::pstderr;
        redi::ipstream child(cmd, mode);
        char buf[1024];
        std::streamsize n;
        bool finished[2] = { false, false };
        while (!finished[0] || !finished[1]) {
            if (!finished[0]) {
                while ((n = child.err().readsome(buf, sizeof(buf))) > 0)
                    err_file.write(buf, n);
                if (child.eof()) {
                    finished[0] = true;
                    if (!finished[1])
                        child.clear();
                }
            }
            if (!finished[1]) {
                while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                    out_file.write(buf, n).flush();
                if (child.eof()) {
                    finished[1] = true;
                    if (!finished[0])
                        child.clear();
                }
            }
        }
        child.close();
        out_file.close();
        err_file.close();
        if (child.rdbuf()->exited())
            return child.rdbuf()->status();
        return 1;
    }
    // todo, may want to handle differently
    // TODO change to sending map of flags as command
    int execute_cmd(std::string cmd) {
        const redi::pstreams::pmode mode = redi::pstreams::pstdout|redi::pstreams::pstderr;
        redi::ipstream child(cmd, mode);
        char buf[1024];
        std::streamsize n;
        bool finished[2] = { false, false };
        while (!finished[0] || !finished[1]) {
            if (!finished[0]) {
                while ((n = child.err().readsome(buf, sizeof(buf))) > 0)
                    continue;
                if (child.eof()) {
                    finished[0] = true;
                    if (!finished[1])
                        child.clear();
                }
            }
            if (!finished[1]) {
                while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                    continue;
                if (child.eof()) {
                    finished[1] = true;
                    if (!finished[0])
                        child.clear();
                }
            }
        }
        child.close();
        if (child.rdbuf()->exited())
            return child.rdbuf()->status();
        return 1;
    }

    void verify_state() {
        // check current state, move to next state
    }

    std::string download_file(std::string flag, std::string &path, std::string temp) {
        // TODO libcurl

        std::string ftp_address;
        std::string output_path = path +flag + ".gz";

        if (flag == ENTAP_CONFIG::INPUT_UNIPROT_SWISS) {
            ftp_address = ENTAP_CONFIG::UNIPROT_FTP_SWISS;

        } else {
            throw ExceptionHandler("Invalid uniprot flag", ENTAP_ERR::E_INPUT_PARSE);
        }

        std::string download_command = "wget -O "+ output_path + " " + ftp_address;
        print_msg("Downloading uniprot: " + flag + " database from " +
                  ftp_address + "...");
        int status = execute_cmd(download_command);
        if (status != 0) {
            throw ExceptionHandler("Error in downloading uniprot database", ENTAP_ERR::E_INIT_TAX_DOWN);
        }
        print_msg("File successfully downloaded to: " + output_path);
        return output_path;
    }

    std::string download_file(const std::string &ftp, std::string &out_path) {
        boostFS::path path(out_path);
        std::string download_command = "wget -O "+ out_path + " " + ftp;
        print_msg("Downloading through wget: file from " + ftp + "...");
        execute_cmd(download_command);
        print_msg("Success! File printed to: " + out_path);
        return out_path;
    }

    void decompress_file(std::string file_path, std::string out_path,short flag) {
        print_msg("Decompressing file at: " + file_path);
        std::string unzip_command;
        if (flag == 0){
            unzip_command = "gzip -d " + file_path;
        } else {
            unzip_command = "tar -xzf " + file_path + " -C " + out_path;
        }
        std::string std_out = out_path + "_std";
        int status = execute_cmd(unzip_command,std_out);
        if (status != 0) {
            throw ExceptionHandler("Error in unzipping database at " +
                    file_path, ENTAP_ERR::E_INIT_TAX_DOWN);
        }
        print_msg("File at: " + file_path + " successfully decompressed");
    }

    int update_database(std::string file_path) {
        return 0;
    }

    int get_supported_threads(boost::program_options::variables_map &user_map) {
        unsigned int supported_threads = std::thread::hardware_concurrency();
        int threads;
        if (user_map["threads"].as<int>() > supported_threads) {
            entapInit::print_msg("Specified thread number is larger than available threads,"
                                         "setting threads to " + std::to_string(supported_threads));
            threads = supported_threads;
        } else {
            threads = user_map["threads"].as<int>();
        }
        return threads;
    }

    std::string generate_command(std::unordered_map<std::string,std::string> &map,
        std::string exe_path) {
        std::stringstream ss;
        ss << exe_path << " ";
        for (auto &pair : map)ss << pair.first << " " << pair.second << " ";
        std::string out = ss.str();
        return out;
    }

}
