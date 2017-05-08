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
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/program_options/variables_map.hpp>
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

    void init_entap(boost::program_options::variables_map user_map, std::string exe_path) {

        // todo print map
//        print_input(input_map);

        std::string outpath = current_path.string() + user_map["tag"].as<std::string>();
        boostFS::path bin_dir(exe_path + ENTAP_CONFIG::BIN_PATH);
        boostFS::path data_dir(exe_path + "/databases");
        boostFS::path out_dir(outpath);
        bool out_dir_state = (boostFS::create_directories(out_dir));
        bool bin_dir_state = (boostFS::create_directories(bin_dir));
        bool data_dir_state = (boostFS::create_directories(data_dir));

        // while state != EXIT_STATE
        try {
            state = INIT_TAX;
            init_taxonomic(exe_path);
//            init_uniprot(input_map["U"]);
//            init_diamond_index(input_map["U"],input_map["N"], input_map["d"]);

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

    void init_uniprot(std::string d) {
        // TODO setup go term/interpro... integration, date tag, use bool input
        std::string ftp_address;
        if (d.compare(ENTAP_CONFIG::INPUT_UNIPROT_SWISS) == 0) {
            ftp_address = ENTAP_CONFIG::UNIPROT_FTP_SWISS;
        } else {
            std::string msg_err = "Invalid input";
            throw ExceptionHandler(msg_err, ENTAP_ERR::E_INPUT_PARSE);
        }

        std::string file_path = ENTAP_CONFIG::UNIPROT_BASE_PATH + d + ".fasta.gz";
        std::string uniprot_command = "wget -O "+ file_path + " " + ftp_address;
        print_msg("Downloading uniprot: " + d + " database from " +
            ftp_address + "...");
        std::cout<<uniprot_command <<std::endl;
        redi::ipstream in(uniprot_command);
        in.close();
        int status = in.rdbuf()->status();
        if (status != 0) {
            throw ExceptionHandler("Error in downloading uniprot database", ENTAP_ERR::E_INIT_TAX_DOWN);
        }

        std::string unzip_command = "gzip -d " + file_path;
        redi::ipstream unz(unzip_command);
        unz.close();
        status = unz.rdbuf()->status();
        if (status != 0) {
            throw ExceptionHandler("Error in unzipping Uniprot Database", ENTAP_ERR::E_INIT_TAX_DOWN);
        }
    }

    void init_ncbi(std::string d) {

    }

    void init_database_parse(std::string path) {

    }

    void init_diamond_index(std::string uniprot, std::string ncbi, std::string database) {
        print_msg("Preparing to index database(s) with Diamond...");
        remove(ENTAP_CONFIG::DIAMOND_INDX_OUT_PATH.c_str());
        if (uniprot.compare(ENTAP_CONFIG::INPUT_UNIPROT_NULL) != 0)  {
//            std::string uniprot_path = ENTAP_CONFIG::UNIPROT_BASE_PATH + uniprot + ".fasta";
            std::string uniprot_path = "databases/database_sequences.fasta";
            if (file_exists(uniprot_path)) {
                print_msg("Uniprot Database - " + uniprot + " found.");
                std::string uniprot_index_command = ENTAP_CONFIG::DIAMOND_PATH_EXE + " makedb --in " +
                        uniprot_path + " -d " + "bin/" + "uniprot_"+uniprot;
                print_msg("Beginning to index database with Diamond...");
                if (execute_cmd(uniprot_index_command, ENTAP_CONFIG::DIAMOND_INDX_OUT_PATH)!=0) {
                    std::cerr<<"Error!!" << std::endl;
                }
                print_msg("Success!");
            } else {
                throw ExceptionHandler("Uniprot database file not found", ENTAP_ERR::E_INIT_INDX_DATA_NOT_FOUND);
            }
        }
        // TODO ncbi and random database array
    }

    void print_input(std::unordered_map<std::string, std::string> map) {
        for ( auto it = map.begin(); it != map.end(); ++it )
            std::cout << " " << it->first << ":" << it->second;
        std::cout << std::endl;
    }

    void print_msg(std::string msg) {
        time_t rawtime;
        time(&rawtime);
        std::string date_time = ctime(&rawtime);
        std::ofstream log_file("debug.txt", std::ios::out | std::ios::app);
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
                    std::cerr.write(buf, n);
                if (child.eof()) {
                    finished[0] = true;
                    if (!finished[1])
                        child.clear();
                }
            }
            if (!finished[1]) {
                while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                    std::cout.write(buf, n).flush();
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
}
