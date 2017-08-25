
//*********************** Includes *****************************
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <boost/filesystem/operations.hpp>
#include <chrono>
#include "EntapConfig.h"
#include "ExceptionHandler.h"
#include "boost/program_options.hpp"
#include "EntapGlobals.h"
#include "EntapExecute.h"
#include "UserInput.h"

//**************************************************************


//******************** Prototype Functions *********************
void init_log();
void init_entap(boostPO::variables_map&);

//**************************************************************


enum States {
    PARSE_ARGS            = 0x01,
    CONFIG_ENTAP          = 0x02,
    CONFIG_ENTAP_SUCCESS  = 0x04,
    EXECUTE_ENTAP         = 0x08,
    EXECUTE_ENTAP_SUCCESS = 0x16
};


States state;   // might handle state summary
std::string _outpath;
std::string _exe_path;
std::string _working_dir;
std::chrono::time_point<std::chrono::system_clock> _start_time;
std::chrono::time_point<std::chrono::system_clock> _end_time;
std::string DEBUG_FILE_PATH;
std::string LOG_FILE_PATH;

int main(int argc, const char** argv) {

    boostPO::variables_map                      inputs;
    std::pair<bool,boostPO::variables_map>      user_pair;

    _start_time = std::chrono::system_clock::now();
    try {
        user_pair = entap_user_parse(argc,argv);
        inputs = user_pair.second;
        init_entap(inputs);
        if (user_pair.first) {
            state = CONFIG_ENTAP;
            entapConfig::init_entap(inputs, _exe_path);
            state = CONFIG_ENTAP_SUCCESS;
        } else {
            state = EXECUTE_ENTAP;
            entapExecute::execute_main(inputs);
            state = EXECUTE_ENTAP_SUCCESS;
        }
    } catch (ExceptionHandler &e) {
        if (e.getErr_code()==ENTAP_ERR::E_SUCCESS) return 0;
        e.print_msg();
        return e.getErr_code();
    }
    _end_time = std::chrono::system_clock::now();
    return 0;
}


void init_entap(boostPO::variables_map& user_input) {

    std::string config_path;

    boost::filesystem::path working_dir(boost::filesystem::current_path());
    _working_dir = working_dir.string();
    _outpath = PATHS(_working_dir, user_input["tag"].as<std::string>());
    DEBUG_FILE_PATH = PATHS(_outpath, ENTAP_CONFIG::DEBUG_FILENAME);
    LOG_FILE_PATH   = PATHS(_outpath, ENTAP_CONFIG::LOG_FILENAME);
    init_log();
    _exe_path = get_exe_path(user_input);
    print_user_input(user_input, _exe_path, _outpath);
    if (user_input.count(ENTAP_CONFIG::INPUT_FLAG_EXE_PATH)) {
        config_path = user_input[ENTAP_CONFIG::INPUT_FLAG_EXE_PATH].as<std::string>();
    } else config_path = PATHS(_exe_path,CONFIG_FILE);
    parse_config(config_path,_exe_path);
}


void init_log() {
    boostFS::remove(DEBUG_FILE_PATH);
    boostFS::remove(LOG_FILE_PATH);
    print_debug("Start - EnTAP");
}
