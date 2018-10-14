/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
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
#include "common.h"
#include <boost/filesystem/operations.hpp>
#include <chrono>
#include <boost/date_time/time_clock.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "EntapConfig.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
#include "EntapExecute.h"
#include "UserInput.h"
#include "EntapConfig.h"
//**************************************************************


//******************** Local Variables *************************
std::chrono::time_point<std::chrono::system_clock> _start_time;
std::chrono::time_point<std::chrono::system_clock> _end_time;
std::string DEBUG_FILE_PATH;        // Extern
std::string LOG_FILE_PATH;          // Extern
UserInput *_pUserInput;
FileSystem *_pFileSystem;
bool _is_config;
//**************************************************************

//******************** Prototype Functions *********************
void init_entap(int, const char**);
void exit_print();

//**************************************************************


/**
 * ======================================================================
 * Function int main(int argc, const char** argv)
 *
 * Description          - Typical C entry point
 *                      - Initializes logs, inputs
 *                      - Sets EnTAP to execution or configuration
 *
 * Notes                - None
 *
 * @param argc          - User input size
 * @param argv          - User input
 * @return              - Execution status, 0 if success
 * ======================================================================
 */
int main(int argc, const char** argv) {
    try {
        init_entap(argc, argv);     // set up logging/user input
        if (_is_config) {
            entapConfig::execute_main(_pUserInput, _pFileSystem);
        } else {
            entapExecute::execute_main(_pUserInput, _pFileSystem);
        }
        exit_print();
    } catch (ExceptionHandler &e) {
        if (e.getErr_code()==ERR_ENTAP_SUCCESS) return 0;
        e.print_msg(_pFileSystem);
        exit_print();   // Might do something different later...
        return e.getErr_code();
    }
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
void init_entap(int argc, const char** argv) {

    pair_str_t  config_default;     // first:config file path, second: default exes
    std::string root_outfiles;
    std::string current_dir;

    // Begin timing
    _start_time = std::chrono::system_clock::now();

    // parse user flags and turn into map
    _pUserInput = new UserInput(argc, argv);

    root_outfiles = _pUserInput->get_user_input<std::string>(_pUserInput->INPUT_FLAG_TAG);

    // create filesystem and begin logging
    _pFileSystem = new FileSystem(root_outfiles);
    _pUserInput->set_pFileSystem(_pFileSystem);

    // get config file path and default to find exe paths
    // First: config path, Second: Software directory defaults
    config_default = _pUserInput->get_config_path();

    // Parse config file
    _pUserInput->parse_config(config_default);

    // Verify and print user input, sets if user selected config or execute
    _is_config = _pUserInput->verify_user_input();
}



/**
 * ======================================================================
 * Function void exit_print(States s)
 *
 * Description          - Prints final information on EnTAP exit
 *
 * Notes                - Will be used for state debug reports in future
 *
 * @return              - None
 * ======================================================================
 */
void exit_print() {
    std::stringstream out_stream;
    std::string       out_msg;
    std::string       out_time;
    long              min_dif;

    _end_time = std::chrono::system_clock::now();
    FS_dprint("End - EnTAP");
    min_dif = std::chrono::duration_cast<std::chrono::minutes>(_end_time - _start_time).count();

    out_stream <<
               "\nEnTAP has completed! "           <<
               "\nTotal runtime (minutes): "       << min_dif;
    out_msg = out_stream.str();
    _pFileSystem->print_stats(out_msg);

    // Cleanup
    delete _pFileSystem;
    delete _pUserInput;
}