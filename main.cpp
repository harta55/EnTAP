#include <iostream>
#include <fstream>
#include <array>
#include <cstring>
#include "EntapInit.h"
#include "ExceptionHandler.h"
#include "pstream.h"
#include "boost/program_options.hpp"
#include "ErrorFlags.h"

namespace boostPO = boost::program_options;

void set_state(int);
void parse_arguments(char**, int);
void print_help();
void print_msg(std::string, bool);
void init_log();
void parse_databases(std::string);
void init_uniprot();
void parse_arguments_boost(int,const char**);

namespace States {
    const static int PARSE_ARGS = 1;
    const static int PARSE_ARGS_SUCCESS = 2;
}

enum STATES_CONFIG {
    PARSE_ARGS          = 0x01
};

STATES_CONFIG state = PARSE_ARGS;   // init

int main(int argc, const char** argv) {
    init_log();

    try {
        set_state(States::PARSE_ARGS);
        parse_arguments_boost(argc,argv);
        set_state(States::PARSE_ARGS_SUCCESS);
//        entapInit::init_entap();
    } catch (ExceptionHandler &e) {
        if (e.getErr_code()==ENTAPERR::E_SUCCESS) return 0;
        e.print_msg();
        return 1;
    }
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
                ("help,h", "help options")
                ("config", "Configure enTAP for execution later")
                ("run", "Execute enTAP functionality")
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
        posOptions.add("run", 1);
        boostPO::variables_map vm;

        try {
            boostPO::store(boostPO::command_line_parser(argc,argv).options(description)
                .positional(posOptions).run(),vm);
            boostPO::notify(vm);

            if (vm.count("help")) {
                std::cout << description<<std::endl<<std::endl;
                throw(ExceptionHandler("",ENTAPERR::E_SUCCESS));
            }

            bool is_config = (bool) vm.count("config");     // ignore 'config config'
            bool is_run = (bool) vm.count("run");

            if (!is_config && !is_run) {
                std::string msg = "Either config option or run option are required";
                throw(ExceptionHandler(msg.c_str(),ENTAPERR::E_INPUT_PARSE));
            }

            if (is_config && is_run) {
                std::string msg = "Cannot specify both config and run flags";
                throw(ExceptionHandler(msg.c_str(),ENTAPERR::E_INPUT_PARSE));
            }
        } catch (boost::program_options::required_option& e) {
            std::cout<<"Required Option"<<std::endl;
        }
    }catch (boost::program_options::error& e){
        // Unknown input
        throw ExceptionHandler(e.what(),ENTAPERR::E_INPUT_PARSE);
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

void init_log() {
    remove("debug.txt");
    print_msg("Start - enTAP", true);
}

// true for additional error pipe, false otherwise
void print_msg(std::string msg, bool b) {
    time_t rawtime;
    time(&rawtime);
    std::string date_time = ctime(&rawtime);
    std::ofstream log_file("debug.txt", std::ios::out | std::ios::app);
    log_file << date_time.substr(0, date_time.size() - 2)
                 + ": " + msg << std::endl;
    log_file.close();
}

void set_state(int flag) {
    state = flag;
}

void set_state(int flag, int next) {
    s
}