/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
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
#include "GraphingManager.h"
#include "FileSystem.h"
#include "TerminalCommands.h"
//**************************************************************

/**
 * ======================================================================
 * Function GraphingManager::GraphingManager(std::string path)
 *
 * Description          - Constructor for graphing manager
 *                      - Checks whether graphing is supported on system
 *
 * Notes                - None
 *
 * @param path          - Path to python graphing file (in /src)
 *
 * @return              - GraphingManager object
 *
 * =====================================================================
 */
GraphingManager::GraphingManager(std::string path) {
    TerminalData terminalData;
    std::string cmd;            // Test command

    FS_dprint("Spawn Object - GraphingManager");

    _graph_path = path;
    cmd = "python " + path + " -s -1 -g -1 -i /temp -t temp";

    terminalData.command = cmd;
    terminalData.print_files = false;

    _graphing_enabled = TC_execute_cmd(terminalData) == 0;
    if (_graphing_enabled) {
        FS_dprint("Graphing is supported");
    } else FS_dprint("Graphing is NOT supported");
}


/**
 * ======================================================================
 * Function void GraphingManager::graph(GraphingStruct& graphingStruct)
 *
 * Description          - Responsible for sending a command to the python
 *                        graphing script
 *
 * Notes                - None
 *
 * @param graphingStruct- Structure of graphing commands
 *
 * @return              - None
 *
 * =====================================================================
 */
void GraphingManager::graph(GraphingData& graphingStruct) {
    if (!_graphing_enabled) return;

    std::unordered_map<std::string,std::string>     cmd_map;
    std::string                                     graphing_cmd;
    TerminalData                                    terminalData;

    cmd_map[FLAG_GRAPH_TEXT] = graphingStruct.text_file_path;
    cmd_map[FLAG_TITLE]      = graphingStruct.graph_title;
    cmd_map[FLAG_OUT_PATH]   = graphingStruct.fig_out_path;
    cmd_map[FLAG_SOFT]       = std::to_string(graphingStruct.software_flag);
    cmd_map[FLAG_GRAPH]      = std::to_string(graphingStruct.graph_type);

    graphing_cmd = generate_command(cmd_map, "python " + _graph_path);

    terminalData.command = graphing_cmd;
    terminalData.print_files = false;

    if (TC_execute_cmd(terminalData) != 0) {
        FS_dprint("\nError generating graph from:\n" + graphing_cmd);
    }
}


bool GraphingManager::is_graphing_enabled() const {
    return _graphing_enabled;
}
