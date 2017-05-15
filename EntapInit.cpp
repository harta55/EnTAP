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
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/program_options/variables_map.hpp>
#include <thread>
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"

//#include "boost/iostreams/filter/gzip.hpp"
//#include "boost/iostreams/copy.hpp"
//#include "boost/iostreams/filtering_streambuf.hpp"
//#include <curl/curl.h> TODO Integrate curllib

namespace boostFS = boost::filesystem;
namespace boostAR = boost::archive;
//namespace boostIO = boost::iostreams;

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

    // may handle differently than ncbi with formatting
    void init_uniprot(std::vector<std::string> &flags, std::string exe) {
        // TODO setup go term/interpro... integration, date tag, use bool input
        print_msg("Parsing uniprot databases...");
        if (flags.empty()) return;
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
                    decompress_file(temp_path);
                    _compiled_databases.push_back(database_path);
                } catch (ExceptionHandler &e) {throw e;}
            }
        }
    }

    void init_ncbi(std::vector<std::string> &flags, std::string exe) {
        // TODO setup go term/interpro... integration, date tag, use bool input
        print_msg("Parsing NCBI databases...");
        if (flags.empty()) return;
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
                    decompress_file(temp_path);
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
        time_t rawtime;
        time(&rawtime);
        std::string date_time = ctime(&rawtime);
        std::ofstream log_file(ENTAP_CONFIG::DEBUG_FILENAME, std::ios::out | std::ios::app);
        log_file << date_time.substr(0, date_time.size() - 2)
                    + ": " + msg << std::endl;
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
                    // std::cerr.write(buf, n);
                if (child.eof()) {
                    finished[0] = true;
                    if (!finished[1])
                        child.clear();
                }
            }
            if (!finished[1]) {
                while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                    continue;
                    // std::cout.write(buf, n).flush();
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

    std::string download_file(std::string flag, std::string &path) {
        std::string ftp_address;
        std::string output_path;

        if (flag == ENTAP_CONFIG::INPUT_UNIPROT_SWISS) {
            ftp_address = ENTAP_CONFIG::UNIPROT_FTP_SWISS;
            output_path += ".gz";

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

    void decompress_file(std::string file_path) {
        int status;
        std::string unzip_command = "gzip -d " + file_path;
        status = execute_cmd(unzip_command);
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
}
