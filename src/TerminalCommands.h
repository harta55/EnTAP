/*
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
*/

#ifndef ENTAP_TERMINALCOMMANDS_H
#define ENTAP_TERMINALCOMMANDS_H

#include "common.h"

#define TC_NULL_ARGUMENT ""

struct TerminalData{
    std::string command;
    std::string out_stream;
    std::string err_stream;
    bool print_files;
    bool suppress_std_err;  // Suppress std error to the debug file
    std::string base_std_path;
};

typedef enum {
    TC_PRINT_COUT,
    TC_PRINT_CERR
} TC_PRINT_TYPES;

typedef std::unordered_map<std::string, std::string> command_map_t;

const int TC_EXIT_SUCCESS = EXIT_SUCCESS;

int TC_execute_cmd(TerminalData &terminalData);
void TC_print(TC_PRINT_TYPES printTypes, std::string msg);
std::string TC_generate_command(command_map_t& command_map, std::string& exe_path);

#endif //ENTAP_TERMINALCOMMANDS_H
