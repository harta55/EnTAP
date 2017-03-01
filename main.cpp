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

enum states_config {
    PARSE_ARGS          = 0x01,
    PARSE_ARGS_SUCCESS  = 0x02,
    INIT_ENTAP          = 0x04,
    INIT_ENTAP_SUCCESS  = 0x08
};

void parse_arguments(char**, int);
void print_help();
void print_msg(std::string, bool);
void init_log();
void parse_arguments_boost(int,const char**);
void state_summary(states_config);

int main(int argc, const char** argv) {
    init_log();
    states_config state_config;   // init

    try {
        state_config = PARSE_ARGS;
        parse_arguments_boost(argc,argv);
//        entapInit::init_entap();
    } catch (ExceptionHandler &e) {
        if (e.getErr_code()==ENTAPERR::E_SUCCESS) return 0;
        e.print_msg();
        state_summary(state_config);
        return 1;
    }
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

void state_summary(states_config st) {
    switch (st) {
        case(PARSE_ARGS):
            break;
        default:
            break;
    }
}