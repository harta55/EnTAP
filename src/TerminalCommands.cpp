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

#include <ios>
#include <pstream.h>
#include "TerminalCommands.h"
#include "common.h"
#include "FileSystem.h"

/**
 * ======================================================================
 * Function int execute_cmd(std::string cmd, std::string out_path)
 *
 * Description          - Terminal stream based on pstreams implementation
 *                      - Executes commands and prints .err and .out from stream
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
int TC_execute_cmd(std::string cmd, std::string out_path) {
    FS_dprint("Executing command: \n" + cmd);
    std::ofstream out_file(out_path+".out", std::ios::out | std::ios::app);
    std::ofstream err_file(out_path+".err", std::ios::out | std::ios::app);
    const redi::pstreams::pmode mode = redi::pstreams::pstdout|redi::pstreams::pstderr;
    redi::ipstream child(cmd, mode);
    char buf[1024];
    std::streamsize n;
    bool finished[2] = { false, false };
    while (!finished[0] || !finished[1]) {
        if (!finished[0]) {
            while ((n = child.err().readsome(buf, sizeof(buf))) > 0)
                err_file.write(buf, n);
            if (child.eof()) {
                finished[0] = true;
                if (!finished[1])
                    child.clear();
            }
        }
        if (!finished[1]) {
            while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                out_file.write(buf, n).flush();
            if (child.eof()) {
                finished[1] = true;
                if (!finished[0])
                    child.clear();
            }
        }
    }
    child.close();
    out_file.close();
    err_file.close();
    if (child.rdbuf()->exited())
        return child.rdbuf()->status();
    return 1;
}
// todo, may want to handle differently
// TODO change to sending map of flags as command
int TC_execute_cmd(std::string cmd) {
    FS_dprint("Executing command: \n" + cmd);
    const redi::pstreams::pmode mode = redi::pstreams::pstdout|redi::pstreams::pstderr;
    redi::ipstream child(cmd, mode);
    char buf[1024];
    std::streamsize n;
    bool finished[2] = { false, false };
    while (!finished[0] || !finished[1]) {
        if (!finished[0]) {
            while ((n = child.err().readsome(buf, sizeof(buf))) > 0)
                continue;
            if (child.eof()) {
                finished[0] = true;
                if (!finished[1])
                    child.clear();
            }
        }
        if (!finished[1]) {
            while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                continue;
            if (child.eof()) {
                finished[1] = true;
                if (!finished[0])
                    child.clear();
            }
        }
    }
    child.close();
    if (child.rdbuf()->exited())
        return child.rdbuf()->status();
    return 1;
}
