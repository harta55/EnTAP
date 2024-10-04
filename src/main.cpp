/******************************************************************
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2024, Alexander Hart, Dr. Jill Wegrzyn
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
 *******************************************************************/


//*********************** Includes *****************************
#include "common.h"
#include <chrono>
#include "config.h"

#ifdef USE_BOOST
#include <boost/filesystem/operations.hpp>
#include <boost/date_time/time_clock.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#endif

#include "EntapConfig.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
#include "EntapExecute.h"
#include "UserInput.h"
#include "EntapConfig.h"
#ifdef UNIT_TESTS
#include "tests/UnitTests.h"
#endif
//********************************************1******************


//******************** Local Variables *************************
UserInput *pUserInput;                                              // Pointer to user input flags
FileSystem *pFileSystem;                                            // Pointer to EnTAP filesystem
UserInput::EXECUTION_TYPE EXECUTE_TYPE;                             // Config/Execute/API
//**************************************************************

//******************** Global Variables ************************
std::string DEBUG_FILE_PATH;        // Extern
std::string LOG_FILE_PATH;          // Extern
bool INITIALIZED_DEBUG_FILE;        // Extern
//**************************************************************

//******************** Local Prototype Functions ***************
void init_entap(int argc, const char** argv);
void exit_print(bool in_error, UserInput::EXECUTION_TYPE);
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
#ifdef UNIT_TESTS
        UnitTests unitTests = UnitTests();
        unitTests.execute_tests();
#else
        init_entap(argc, argv);     // set up logging/user input
        if (EXECUTE_TYPE == UserInput::EXECUTE_CONFIG) {
            entapConfig::execute_main(pUserInput, pFileSystem);
        } else if (EXECUTE_TYPE == UserInput::EXECUTE_EXECUTE) {
            entapExecute::execute_main(pUserInput, pFileSystem);
        } else if (EXECUTE_TYPE == UserInput::EXECUTE_API) {
            ;
        } else {
            ;
        }
        exit_print(false, EXECUTE_TYPE);
        SAFE_DELETE(pFileSystem);
        SAFE_DELETE(pUserInput);
#endif

    } catch (ExceptionHandler &e) {
        if (e.getErr_code()==ERR_ENTAP_SUCCESS) return 0;
        e.print_msg(pFileSystem);
        exit_print(true, EXECUTE_TYPE);
        return e.getErr_code();
    }
    return 0;
}


/**
 * ======================================================================
 * Function init_entap(int argc, const char** argv)
 *
 * Description          - Initializes debug/log file paths
 *                      - Parses user configuration file and inputs
 *                      - Initializes file system
 *
 * Notes                - None
 *
 * @param argc          - User input size
 * @param argv          - User input
 *
 * @return              - None
 * ======================================================================
 */
void init_entap(int argc, const char** argv) {

    // Begin timing
    auto const startTime = std::chrono::system_clock::now();
    TC_print(TC_PRINT_COUT, get_cur_time() + " -- Start EnTAP");
    INITIALIZED_DEBUG_FILE = false;

    // Create filesystem
    pFileSystem = new FileSystem(startTime);

    // parse user flags and set the root filesystem directory according to user input
    pUserInput = new UserInput(argc, argv, pFileSystem);

    // Verify and print user input (marks begining of log file), sets if user selected config or execute
    EXECUTE_TYPE = pUserInput->verify_user_input();
}


/**
 * ======================================================================
 * Function void exit_print(bool in_error)
 *
 * Description          - Prints final information on EnTAP exit
 *
 * Notes                - Kills FileSystem and UserInput objects
 *
 * @param in_error      - Indicates whether the system is finishing in error
 *                        (true) or not (false)
 *
 * @return              - None
 * ======================================================================
 */
void exit_print(bool in_error, UserInput::EXECUTION_TYPE execute_type) {

    switch (execute_type) {
        // FALL THROUGH for EXECUTE_EXECUTE
        case UserInput::EXECUTE_EXECUTE:
        case UserInput::EXECUTE_CONFIG: {
            std::stringstream out_stream;  // Stream to print to EnTAP log file
            std::string out_msg;     // Temp string of outstream
            int64 min_dif;     // start time to end time difference
            std::chrono::time_point<std::chrono::system_clock> end_time;    // EnTAP end time

            end_time = std::chrono::system_clock::now();
            FS_dprint("End - EnTAP");
            min_dif = std::chrono::duration_cast<std::chrono::minutes>(end_time - pFileSystem->m_start_time()).count();

            if (in_error) {
                out_stream <<
                           "\nEnTAP has failed! ";
            } else {
                out_stream <<
                           "\nEnTAP has completed!";
            }
            out_stream << "\nTotal runtime (minutes): " << min_dif;
            out_msg = out_stream.str();
            pFileSystem->print_stats(out_msg);
            TC_print(TC_PRINT_COUT, get_cur_time() + " -- EnTAP Complete! [total runtime " + std::to_string(min_dif) +
                " minute(s)]");
            break;
        }

        case UserInput::EXECUTE_API:
        {
            std::string out_msg;
            out_msg = pUserInput->get_json_output();
            TC_print(TC_PRINT_COUT, out_msg);
            break;
        }

        default:
            break;
    }
}