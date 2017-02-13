//
// Created by harta55 on 2/1/17.
//

#include <cstdio>
#include "EntapInit.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include "pstream.h"
#include "boost/filesystem.hpp"
//#include "boost/iostreams/filter/gzip.hpp"
//#include "boost/iostreams/copy.hpp"
//#include "boost/iostreams/filtering_streambuf.hpp"
//#include <curl/curl.h> TODO Integrate curllib

namespace boostFS = boost::filesystem;
//namespace boostIO = boost::iostreams;

int init_taxonomic();

namespace entapInit {

    const boostFS::path current_path (boost::filesystem::current_path());


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


    void init_entap() {
        std::cout<<"hello"<<std::endl;

        boostFS::path bin_dir(current_path.string()+"/bin");
        boostFS::path data_dir(current_path.string()+"/databases");
        bool bin_dir_state = (boostFS::create_directories(bin_dir));
        bool data_dir_state = (boostFS::create_directories(data_dir));
        init_taxonomic();



    }

    bool file_exists (const std::string& name) {
        struct stat buff;
        return (stat (name.c_str(), &buff) == 0);
    }

    int init_taxonomic() {
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
        std::string database_path = current_path.string()+"/databases";
        std::vector<std::string> errors;
        redi::ipstream in("wget ftp://ftp.pir.georgetown.edu/databases/idmapping/refseq/uniprotkb2entrez_gene.gz "
                                  "-P "+ database_path,
                          redi::pstreambuf::pstderr);
        std::string errmsg;
        std::string str;
        while (std::getline(in>>str, errmsg)) {
            errors.push_back(errmsg);
            std::cout<<errmsg<<std::endl;
            std::cout<<str<<std::endl;
        }
        if (errors.empty()) std::cout<<"NO ERROR"<<std::endl;

    }

//    int decompress_gzip()
}
