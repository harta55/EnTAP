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
#include "ErrorFlags.h"
#include "ExceptionHandler.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_set.hpp>
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"

//#include "boost/iostreams/filter/gzip.hpp"
//#include "boost/iostreams/copy.hpp"
//#include "boost/iostreams/filtering_streambuf.hpp"
//#include <curl/curl.h> TODO Integrate curllib

namespace boostFS = boost::filesystem;
namespace boostAR = boost::archive;
//namespace boostIO = boost::iostreams;

int init_taxonomic();

namespace entapInit {

    const boostFS::path current_path(boost::filesystem::current_path());


//    size_t write_function(void *buffer, size_t size, size_t nmemb, void *stream) {
//        struct FtpFile *out = (struct FtpFile *) stream;
//        if (out && !out -> stream) {
//            // open file to write
//            out -> stream = fopen(out->filename, "wb");
//            if (!out->stream) return 1; // can't open file to write
//        }
//        return fwrite(buffer, size, nmemb, out->stream);
//    }
//
//    void download_file(std::string, std::string) {
//        CURL *curl;
//        CURLcode res;
//        struct FtpFile ftpfile={
//                "file.bin", /* name to store the file as if successful */
//                NULL
//        };
//        curl_global_init(CURL_GLOBAL_DEFAULT);
//        curl = curl_easy_init();
//        if(curl) {
//            curl_easy_setopt(curl, CURLOPT_URL,
//                             "ftp_url");
//            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_function);
//            curl_easy_setopt(curl, CURLOPT_WRITEDATA, &ftpfile);
//            curl_easy_setopt(curl, CURLOPT_USE_SSL, CURLUSESSL_ALL);
//            curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);
//            res = curl_easy_perform(curl);
//            curl_easy_cleanup(curl);
//
//            if(CURLE_OK != res) {
//                /* we failed */
//                fprintf(stderr, "curl told us %d\n", res);
//            }
//        }
//        if(ftpfile.stream)
//            fclose(ftpfile.stream); /* close the local file */
//        curl_global_cleanup();
//    }

    void init_entap(std::unordered_map<std::string, std::string> input_map) {

        boostFS::path bin_dir(current_path.string() + "/bin");
        boostFS::path data_dir(current_path.string() + "/databases");
        bool bin_dir_state = (boostFS::create_directories(bin_dir));
        bool data_dir_state = (boostFS::create_directories(data_dir));

        try {
//            init_taxonomic();
            init_uniprot(input_map["a"]);

        }catch (ExceptionHandler &e) {
            throw ExceptionHandler(e.what(), e.getErr_code());
        }
    }

    bool file_exists(const std::string &name) {
        struct stat buff;
        return (stat(name.c_str(), &buff) == 0);
    }

    void init_taxonomic() {
        print_msg("Downloading taxonomic database...");
        //TODO Integrate gzip/zlib
//        std::ifstream file("file.gz", std::ios_base::in | std::ios_base::binary);
//        try {
//            boostIO::filtering_streambuf<boostIO::input> in;
//            in.push(boostIO::gzip_decompressor());
//            in.push(file);
//
//
//        }
//        catch(const boost::iostreams::gzip_error& e) {
//            std::cout << e.what() << '\n';
//        }

        std::string database_path = current_path.string() + "/databases";
        std::string tax_command = "perl " + ENTAPERR::taxonomic_script;
        redi::ipstream in(tax_command);
        in.close();
        int status = in.rdbuf()->status();
        if (status != 0) {
            std::cerr << "Error in downloading taxonomic database" << std::endl;
            throw ExceptionHandler("Error in downloading taxonomic database", ENTAPERR::E_INIT_TAX_DOWN);
        }
        print_msg("Success!");
        print_msg("Indexing taxonomic database...");

        std::unordered_map<std::string, std::string> tax_data_map;
        std::ifstream infile(ENTAPERR::taxonomic_database);
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
            throw ExceptionHandler(e.what(), ENTAPERR::E_INIT_TAX_INDEX);
        }

        // write map
        try{
            {
                std::ofstream ofs("test.entp");
                boostAR::binary_oarchive oa(ofs);
                oa << tax_data_map;
            }
        } catch (std::exception &e) {
            throw ExceptionHandler(e.what(), ENTAPERR::E_INIT_TAX_SERIAL);
        }

        // read map
//        std::unordered_map<std::string, std::string> restored_map;
//        try {
//            {
//                std::ifstream ifs("test.ent");
//                boost::archive::binary_iarchive ia(ifs);
//                ia >> restored_map;
//            }
//
//            std::cout<<restored_map["pseudoascotaiwania"];
//        } catch (std::exception &exception){
//            throw ExceptionHandler(exception.what(), ENTAPERR::E_INIT_TAX_READ);
//        }
        print_msg("Success!");
    }

    void init_uniprot(std::string d) {
        // TODO setup go term/interpro... integration, date tag
        std::string base = "www.uniprot.org/";
        std::string query;

        if (d.compare("swiss")) {
            std::cout <<  "correct!"<< std::endl;

        }
    }

    void init_ncbi(std::string d) {

    }

    void init_database_parse(std::string path) {

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
}
