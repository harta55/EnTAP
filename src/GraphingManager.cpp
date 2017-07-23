//
// Created by harta on 7/21/17.
//

#include "GraphingManager.h"
#include "EntapInit.h"

GraphingManager::GraphingManager(std::string path) {
    _graph_path = path;
    std::string cmd = "python " + path + " -s -1 -g -1 -i /temp";
    _graphing_enabled = entapInit::execute_cmd(cmd) == 0;
    if (_graphing_enabled) {
        entapInit::print_msg("Graphing is supported");
    } else entapInit::print_msg("Graphing is NOT supported");
}

void GraphingManager::graph(std::string txt_path,std::string title, std::string fig_path,
                            unsigned char soft_flag, unsigned char graph_flag) {
    if (!_graphing_enabled) return;

    std::unordered_map<std::string,std::string>     cmd_map;
    std::string                                     graphing_cmd;

    cmd_map[FLAG_GRAPH_TEXT] = txt_path;
    cmd_map[FLAG_TITLE]      = title;
    cmd_map[FLAG_OUT_PATH]   = fig_path;
    cmd_map[FLAG_SOFT]       = std::to_string(soft_flag);
    cmd_map[FLAG_GRAPH]      = std::to_string(graph_flag);

    graphing_cmd = entapInit::generate_command(cmd_map, "python " + _graph_path);
    if (entapInit::execute_cmd(graphing_cmd) != 0) {
        entapInit::print_msg("\nError generating graph from:\n" + graphing_cmd);
    }
}
