#include <iostream>
#include <fstream>
#include <array>
#include <cstring>
#include "ExceptionHandler.h"
#include "pstream.h"
#include "boost/program_options.hpp"
#include "ErrorFlags.h"

namespace boostPO = boost::program_options;

void parse_arguments(char**, int);
void print_help();
void print_msg(std::string, bool);
void init_log();
void parse_databases(std::string);
void init_uniprot();
void parse_arguments_boost(int,const char**);

int main(int argc, const char** argv) {
//    redi::ipstream in("ls");
//    std::string str;
//    while (in >> str) {
//        std::cout << str << std::endl;
//        std::cout << "database" << std::endl;
//    }
    init_log();
    parse_arguments_boost(argc,argv);
//    try {
//        parse_arguments(argv, argc);
//    }catch (ExceptionHandler e) {
//        switch (e.getErrTag()) {
//            case ExceptionHandler::except_help:
//                return 0;
//            case ExceptionHandler::except_input_parse:
//                print_msg("End - enTAP", true);
//                print_msg("End - Input parse", true);
//                return 1;
//            default:
//                return 1;
//        }
//    }
    return 0;
}

void parse_arguments_boost(int argc, const char** argv) {
    try {
        boostPO::options_description description("Options");
        // TODO separate out into main options and additional with defaults
        description.add_options()
                ("help,h", "Welcome to enTAP!")
                ("config", "Configure enTAP for execution later")
                ("blast", "Execute enTAP functionality")
                ("ncbi,N", boostPO::value<std::string>(),"Select which NCBI database you would like to download"
                        "\nref - RefSeq database...")
                ("uniprot,U", "Select which Uniprot database you would like to download"
                        "\n100 - UniRef100...")
                //multiple entries
                ("database,d", boostPO::value<std::vector<std::string>>(),
                        "Provide the path to a separate database, however this "
                        "may prohibit taxonomic filtering.")
                ("version,v", "Display version number");
        boostPO::positional_options_description posOptions;
        posOptions.add("config", 1);
        posOptions.add("blast", 1);
        boostPO::variables_map vm;

        try {
            boostPO::store(boostPO::command_line_parser(argc,argv).options(description)
                .positional(posOptions).run(),vm);


            if (vm.count("help")) {
                std::cout << description<<std::endl<<std::endl;
            }
            boostPO::notify(vm);
        } catch (boost::program_options::required_option& e) {
            std::cerr <<"ERROR"<<std::endl;
        }
    }catch (ExceptionHandler &e){
        std::cerr <<"ERROR2"<< e.what()<<std::endl;
    }
}

//void parse_arguments(char** args, int len) {
//    print_msg("Start - parse arguments", true);
//    unsigned char input = 0x0;
//    enum Inputs {
//        ARG_CONFIG          = 0x01,
//        ARG_HELP            = 0x02,
//        ARG_OTHER           = 0x04
//    };
//
//    //-c config (download databases / rapsearch)
//    if (len == 1) {
//        // no inputs
//        throw ExceptionHandler(<#initializer#>, ExceptionHandler::except_input_parse, "No inputs");
//    }
//    std::string arg_str = "";
//    for (int i = 1; i != len; i++) {
//        std::string s = args[i];
//        arg_str.append(s + " ");
//        if (s.compare("-h")==0) {
//            input |= ARG_HELP;
//        }else if (s.compare("-c")==0) {
//            input |= ARG_CONFIG;
//        }
//    }
//    if (input & ARG_HELP) {
//        print_help();
//        throw ExceptionHandler(<#initializer#>, ExceptionHandler::except_help, "");
//    } else if (input & ARG_CONFIG & ~ARG_HELP & ~ARG_OTHER) {
//        try {
//            parse_databases(arg_str);
//        }catch (ExceptionHandler e) {
//            throw e;
//        }
//    }
//}

void print_help() {
    // TODO help information
    std::cout<<"Print help"<<std::endl;
}

void init_log() {
    remove("debug.txt");
    print_msg("Start - enTAP", true);
}

// true for additional error pipe, false otherwise
void print_msg(std::string msg, bool b) {
    time_t rawtime;
    time(&rawtime);
    std::string date_time = ctime(&rawtime);
    std::ofstream log_file("debug.txt", std::ios_base::out | std::ios_base::app );
    if (b) {
        log_file << date_time.substr(0, date_time.size() - 2)
                     + ": " + msg << std::endl;
    }
}

//void parse_databases(std::string str) {
//    unsigned char flag = 0x0;
//    enum flags {
//        DATABASE_UNIPROT            = 0x01,
//        DATABASE_NCBI               = 0x02
//    };
//    std::string dbs = str.substr(str.find("-c")+3);
//    if (dbs.length() > 4) {
//        throw ExceptionHandler(<#initializer#>, ExceptionHandler::except_input_parse, "Incorrect config format: " + dbs);
//    }
//    // accept any intermediate character (ie. +,/,-)
//    for (char i:dbs) {
//        if (i == 'U'){
//            flag |= DATABASE_UNIPROT;
//            std::cout << "uniprot" << std::endl;
//        } else if (i == 'N') {
//            flag |= DATABASE_NCBI;
//        }
//    }
//    if (flag & ~(DATABASE_NCBI | DATABASE_UNIPROT)){
//        throw ExceptionHandler(<#initializer#>, ExceptionHandler::except_input_parse, "No databases selected");
//    }
//
//    if (flag & DATABASE_UNIPROT) {
//        init_uniprot();
//    }
//}

void init_uniprot() {
    print_msg("Start - Uniprot download", true);
}