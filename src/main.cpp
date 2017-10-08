/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
 *
 * This file is part of EnTAP.
 *
 * EnTAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EnTAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/


//*********************** Includes *****************************
#include <iostream>
#include <fstream>
#include <array>
#include <ctime>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <boost/filesystem/operations.hpp>
#include <chrono>
#include <iomanip>
#include "EntapConfig.h"
#include "ExceptionHandler.h"
#include "boost/program_options.hpp"
#include "EntapGlobals.h"
#include "EntapExecute.h"
#include "UserInput.h"

//**************************************************************


enum States {
    PARSE_ARGS            = 0x01,
    CONFIG_ENTAP          = 0x02,
    CONFIG_ENTAP_SUCCESS  = 0x04,
    EXECUTE_ENTAP         = 0x08,
    EXECUTE_ENTAP_SUCCESS = 0x16
};

std::string _outpath;
std::string _exe_path;
std::string _working_dir;
std::chrono::time_point<std::chrono::system_clock> _start_time;
std::chrono::time_point<std::chrono::system_clock> _end_time;
std::string DEBUG_FILE_PATH;        // Extern
std::string LOG_FILE_PATH;          // Extern


//******************** Prototype Functions *********************
void init_log();
void init_entap(boostPO::variables_map&);
void exit_print(States);
//**************************************************************


/**
 * ======================================================================
 * Function int main(int argc, const char** argv)
 *
 * Description          - Typical C entry point
 *                      - Initializes logs, inputs
 *                      - Sets EnTAP to execution or configuration states
 *
 * Notes                - None
 *
 * @param argc          - User input size
 * @param argv          - User input
 * @return              - Execution status, 0 if success
 * ======================================================================
 */
int main(int argc, const char** argv) {

    boostPO::variables_map                      inputs;
    States                                      state;
    bool                                        config;

    _start_time = std::chrono::system_clock::now();
    try {
        inputs = parse_arguments_boost(argc,argv);
        init_entap(inputs);
        config = verify_user_input(inputs);
        if (config) {
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
    exit_print(state);
    return 0;
}


/**
 * ======================================================================
 * Function init_entap(boostPO::variables_map& user_input)
 *
 * Description          - Initializes debug/log file paths
 *                      - Parses user configuration file
 *
 * Notes                - None
 *
 * @param user_input    - Boost user inputs
 * @return              - None
 * ======================================================================
 */
void init_entap(boostPO::variables_map& user_input) {

    std::string config_path;

    boost::filesystem::path working_dir(boost::filesystem::current_path());
    _working_dir    = working_dir.string();
    _outpath        = PATHS(_working_dir, user_input["tag"].as<std::string>());
    boostFS::create_directories(_outpath);
    init_log();
    _exe_path = get_exe_path(user_input);
    print_user_input(user_input, _exe_path, _outpath);
    if (user_input.count(ENTAP_CONFIG::INPUT_FLAG_EXE_PATH)) {
        config_path = user_input[ENTAP_CONFIG::INPUT_FLAG_EXE_PATH].as<std::string>();
    } else config_path = PATHS(_exe_path,CONFIG_FILE);
    parse_config(config_path,_exe_path);
}


/**
 * ======================================================================
 * Function init_log()
 *
 * Description          - Initializes log and debug files
 *
 * Notes                - None
 *
 * @return              - None
 * ======================================================================
 */
void init_log() {
    std::chrono::time_point<std::chrono::system_clock> now;
    std::time_t                                        now_time;
    std::tm                                            now_tm;
    std::stringstream                                  ss;
    std::string                                        log_file_name;
    std::string                                        debug_file_name;
    std::string                                        time_date;

    now      = std::chrono::system_clock::now();
    now_time = std::chrono::system_clock::to_time_t(now);
    now_tm   = *std::localtime(&now_time);
    ss << std::put_time(&now_tm, "_%Y.%m.%d-%Hh.%Mm.%Ss");
    time_date       = ss.str();
    log_file_name   = ENTAP_CONFIG::LOG_FILENAME   + time_date + ENTAP_CONFIG::LOG_EXTENSION;
    debug_file_name = ENTAP_CONFIG::DEBUG_FILENAME + time_date + ENTAP_CONFIG::LOG_EXTENSION;
    DEBUG_FILE_PATH = PATHS(_outpath, debug_file_name);
    LOG_FILE_PATH   = PATHS(_outpath, log_file_name);
    boostFS::remove(DEBUG_FILE_PATH);
    boostFS::remove(LOG_FILE_PATH);
    print_debug("Start - EnTAP");
}


/**
 * ======================================================================
 * Function void exit_print(States s)
 *
 * Description          - Prints final information on EnTAP exit
 *
 * Notes                - Will be used for state debug reports in future
 *
 * @param s             - State
 * @return              - None
 * ======================================================================
 */
void exit_print(States s) {
    std::stringstream out_stream;
    std::string       out_msg;
    std::string       out_time;
    long              min_dif;

    print_debug("End - EnTAP");
    min_dif = std::chrono::duration_cast<std::chrono::minutes>(_end_time - _start_time).count();

    out_stream <<
               "\nEnTAP has completed! "           <<
               "\nTotal runtime (minutes): "       << min_dif;
    out_msg = out_stream.str();
    print_statistics(out_msg);
}