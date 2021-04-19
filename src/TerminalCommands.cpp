/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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

#include <pstream.h>
#include "TerminalCommands.h"
#include "common.h"
#include "FileSystem.h"


/**
 * ======================================================================
 * Function int execute_cmd(std::string cmd, std::stringstream err_stream, std::stringstream out_stream)
 *
 * Description          - Terminal stream based on pstreams implementation
 *                      - Executes commands and prints to err/out stream
 *
 * Notes                - None
 *
 * @param cmd           - Command for child process
 * @param out_path      - Path to std out/err files to be printed
 *
 * @return              - int error code
 *
 * =====================================================================
 */
int TC_execute_cmd(TerminalData &terminalData) {

    FS_dprint("Executing command: \n" + terminalData.command);

    std::stringstream err_stream;
    std::stringstream out_stream;

    const redi::pstreams::pmode mode = redi::pstreams::pstdout|redi::pstreams::pstderr;
    redi::ipstream child(terminalData.command, mode);
    char buf[1024];
    std::streamsize n;
    bool finished[2] = { false, false };
    while (!finished[0] || !finished[1]) {
        if (!finished[0]) {
            while ((n = child.err().readsome(buf, sizeof(buf))) > 0)
                err_stream.write(buf, n);
            if (child.eof()) {
                finished[0] = true;
                if (!finished[1])
                    child.clear();
            }
        }
        if (!finished[1]) {
            while ((n = child.out().readsome(buf, sizeof(buf))) > 0) {
                out_stream.write(buf, n).flush();
            }
            if (child.eof()) {
                finished[1] = true;
                if (!finished[0])
                    child.clear();
            }
        }
    }
    child.close();

    terminalData.err_stream = err_stream.str();
    terminalData.out_stream = out_stream.str();

    // Print error to debug file if we are not suppressing output
    if (!terminalData.suppress_std_err) {
        FS_dprint("\nStd Err:\n" + terminalData.err_stream);\
    } else {
        FS_dprint("WARNING: Suppressing error log from child process, check error file generated for full output");
    }

    if (terminalData.print_files) {
        std::string out_path = terminalData.base_std_path + FileSystem::EXT_OUT;
        std::string err_path = terminalData.base_std_path + FileSystem::EXT_ERR;

        std::ofstream out_file(out_path, std::ios::out | std::ios::app);
        std::ofstream err_file(err_path, std::ios::out | std::ios::app);
        FS_dprint("\nPrinting to files:\nStd Out: " + out_path + "\nStd Err: " + err_path);

        out_file << terminalData.out_stream;
        err_file << terminalData.err_stream;

        out_file.close();
        err_file.close();
    }

    if (child.rdbuf()->exited())
        return child.rdbuf()->status();
    return 1;
}

std::string TC_generate_command(command_map_t &map, std::string &exe_path) {
    std::stringstream ss;
    std::string       out;

    ss << exe_path << " ";
    for (auto &pair : map)ss << pair.first << " " << pair.second << " ";
    out = ss.str();
    return out;
}