//
// Created by harta on 7/21/17.
//

#include "GraphingManager.h"
#include "EntapConfig.h"
#include "EntapConsts.h"

GraphingManager::GraphingManager(std::string path) {
    _graph_path = path;
    std::string cmd = "python " + path + " -s -1 -g -1 -i /temp -t temp";
    _graphing_enabled = entapConfig::execute_cmd(cmd) == 0;
    if (_graphing_enabled) {
        entapConfig::print_msg("Graphing is supported");
    } else entapConfig::print_msg("Graphing is NOT supported");
}

void GraphingManager::graph(GraphingStruct& graphingStruct) {
    if (!_graphing_enabled) return;

    std::unordered_map<std::string,std::string>     cmd_map;
    std::string                                     graphing_cmd;

    cmd_map[FLAG_GRAPH_TEXT] = graphingStruct.text_file_path;
    cmd_map[FLAG_TITLE]      = graphingStruct.graph_title;
    cmd_map[FLAG_OUT_PATH]   = graphingStruct.fig_out_path;
    cmd_map[FLAG_SOFT]       = std::to_string(graphingStruct.software_flag);
    cmd_map[FLAG_GRAPH]      = std::to_string(graphingStruct.graph_type);

    graphing_cmd = entapConfig::generate_command(cmd_map, "python " + _graph_path);
    if (entapConfig::execute_cmd(graphing_cmd) != 0) {
        entapConfig::print_msg("\nError generating graph from:\n" + graphing_cmd);
    }
}
