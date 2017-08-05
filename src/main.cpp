
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
std::string get_exe_path(boostPO::variables_map&);

//**************************************************************


enum States {
    PARSE_ARGS            = 0x01,
    INIT_ENTAP            = 0x02,
    INIT_ENTAP_SUCCESS    = 0x04,
    EXECUTE_ENTAP         = 0x08,
    EXECUTE_ENTAP_SUCCESS = 0x16
};


States state;   // init
std::string _outpath;
std::string _exe_path;
std::string _working_dir;
std::chrono::time_point<std::chrono::system_clock> _start_time;
std::chrono::time_point<std::chrono::system_clock> _end_time;
std::string DEBUG_FILE_PATH;
std::string LOG_FILE_PATH;

int main(int argc, const char** argv) {

    std::unordered_map<std::string,std::string> config_map;
    boostPO::variables_map                      inputs;
    std::pair<bool,boostPO::variables_map>      user_pair;
    bool                                        config;     // True for config state

    _start_time = std::chrono::system_clock::now();
    try {
        user_pair = entap_user_parse(argc,argv);
        config = user_pair.first;
        inputs = user_pair.second;
        init_entap(inputs);
        config_map = parse_config(_exe_path);
        print_user_input(inputs, _exe_path, _outpath);
        if (config) {
            state = INIT_ENTAP;
            entapConfig::init_entap(inputs, _exe_path, config_map);
            state = INIT_ENTAP_SUCCESS;
        } else {
            entapExecute::execute_main(inputs, _exe_path, config_map);
            state = EXECUTE_ENTAP_SUCCESS;
        }
    } catch (ExceptionHandler &e) {
        if (e.getErr_code()==ENTAP_ERR::E_SUCCESS) return 0;
        e.print_msg();
        return 1;
    }
    _end_time = std::chrono::system_clock::now();
    return 0;
}


void init_entap(boostPO::variables_map& user_input) {
    boost::filesystem::path working_dir(boost::filesystem::current_path());
    _working_dir = working_dir.string();
    _outpath = (boostFS::path(_working_dir) /
                boostFS::path(user_input["tag"].as<std::string>())).string();
    _exe_path = get_exe_path(user_input);
    DEBUG_FILE_PATH = (boostFS::path(_outpath) / ENTAP_CONFIG::DEBUG_FILENAME).string();
    LOG_FILE_PATH   = (boostFS::path(_outpath) / ENTAP_CONFIG::LOG_FILENAME).string();
    init_log();
}


void init_log() {
    boostFS::remove(DEBUG_FILE_PATH);
    boostFS::remove(LOG_FILE_PATH);
    print_debug("Start - enTAP");
}


std::string get_exe_path(boostPO::variables_map &vm) {
    //TODO check different systems
    if (vm.count(ENTAP_CONFIG::INPUT_FLAG_EXE_PATH)) {
        return vm[ENTAP_CONFIG::INPUT_FLAG_EXE_PATH].as<std::string>();
    }
    char buff[1024];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
        buff[len] = '\0';
        std::string path = std::string(buff);
        boost::filesystem::path p(path);p.remove_filename();
        return p.string();
    }
    return "";
}