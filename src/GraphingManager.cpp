/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/


//*********************** Includes *****************************
#include "GraphingManager.h"
#include "EntapConfig.h"
#include "EntapGlobals.h"
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
    print_debug("Spawn object - GraphingManager");
    _graph_path = path;
    std::string cmd = "python " + path + " -s -1 -g -1 -i /temp -t temp";
    _graphing_enabled = execute_cmd(cmd) == 0;
    if (_graphing_enabled) {
        print_debug("Graphing is supported");
    } else print_debug("Graphing is NOT supported");
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
void GraphingManager::graph(GraphingStruct& graphingStruct) {
    if (!_graphing_enabled) return;

    std::unordered_map<std::string,std::string>     cmd_map;
    std::string                                     graphing_cmd;

    cmd_map[FLAG_GRAPH_TEXT] = graphingStruct.text_file_path;
    cmd_map[FLAG_TITLE]      = graphingStruct.graph_title;
    cmd_map[FLAG_OUT_PATH]   = graphingStruct.fig_out_path;
    cmd_map[FLAG_SOFT]       = std::to_string(graphingStruct.software_flag);
    cmd_map[FLAG_GRAPH]      = std::to_string(graphingStruct.graph_type);

    graphing_cmd = generate_command(cmd_map, "python " + _graph_path);
    if (execute_cmd(graphing_cmd) != 0) {
        print_debug("\nError generating graph from:\n" + graphing_cmd);
    }
}


bool GraphingManager::is_graphing_enabled() const {
    return _graphing_enabled;
}
