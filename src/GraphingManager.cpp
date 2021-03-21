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


//*********************** Includes *****************************
#include "GraphingManager.h"
#include "FileSystem.h"
#include "TerminalCommands.h"
//**************************************************************

/**
 * ======================================================================
 * Function GraphingManager::GraphingManager(std::string path, FileSystem *fileSystem)
 *
 * Description          - Constructor for graphing manager
 *                      - Checks whether graphing is supported on system
 *
 * Notes                - None
 *
 * @param path          - Path to python graphing file (in /src)
 * @param fileSystem    - Pointer to Entap FileSystem object
 *
 * @return              - GraphingManager object
 *
 * =====================================================================
 */
GraphingManager::GraphingManager(std::string path, FileSystem *fileSystem) {
    TerminalData terminalData;
    std::string cmd;            // Test command

    FS_dprint("Spawn Object - GraphingManager");
    mpFileSystem = fileSystem;
    mGraphingPath = path;
    cmd = "python " + path + " " + FLAG_GRAPH_TYPE + " " + std::to_string(ENT_GRAPH_TEST);

    terminalData.command = cmd;
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;

    try {
        mGraphingEnabled = TC_execute_cmd(terminalData) == 0;
        if (mGraphingEnabled) {
            FS_dprint("Graphing is supported");
        } else FS_dprint("Graphing is NOT supported");
    } catch (...) {
        // Disable graphing if any exception thrown
        mGraphingEnabled = false;
        FS_dprint("Graphing is NOT supported, unknown exception");
    }
}

GraphingManager::~GraphingManager() {
    FS_dprint("Killing Object - GraphingManager");
    for (auto &pair : mGraphData) {
        if (pair.second != nullptr) {
            pair.second->close_graphing_file();
            SAFE_DELETE(pair.second);
        }
        mGraphData.erase(pair.first);
    }
}

/**
 * ======================================================================
 * Function void GraphingManager::graph(GraphingStruct* graphingStruct)
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
void GraphingManager::graph(GraphingData* graphingStruct) {
    if (!mGraphingEnabled || graphingStruct == nullptr) return;

    std::unordered_map<std::string,std::string>     cmd_map;
    std::string                                     graphing_cmd;
    std::string                                     exe_path;
    TerminalData                                    terminalData;

    // Python call will change...
    exe_path = "python " + mGraphingPath;
    cmd_map[FLAG_GRAPH_DATA_PATH] = graphingStruct->text_file_path;
    cmd_map[FLAG_TITLE]           = graphingStruct->graph_title;
    cmd_map[FLAG_OUT_PATH]        = graphingStruct->fig_out_path;
    cmd_map[FLAG_GRAPH_TYPE]      = std::to_string(graphingStruct->graph_type);

    graphing_cmd = TC_generate_command(cmd_map, exe_path);

    terminalData.command = graphing_cmd;
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;

    if (TC_execute_cmd(terminalData) != 0) {
        FS_dprint("\nError generating graph from:\n" + graphing_cmd);
    }
}

bool GraphingManager::is_graphing_enabled() const {
    return mGraphingEnabled;
}

bool GraphingManager::initialize_graph_data(GraphingManager::GraphingData &graphingData) {
    FS_dprint("GraphingManager - initializing graphing data at: " + graphingData.text_file_path);
    EntapGraphBase *new_graph = nullptr;
    bool ret = true;

    switch (graphingData.graph_type) {
        case ENT_GRAPH_BAR_VERTICAL:
            new_graph = new EntapGraphBarVertical(graphingData);
            break;

        case ENT_GRAPH_BAR_HORIZONTAL:
            new_graph = new EntapGraphBarHorizontal(graphingData);
            break;

        case ENT_GRAPH_LINE_HORIZONTAL:
            new_graph = new EntapGraphLineHorizontal(graphingData);
            break;

        case ENT_GRAPH_LINE_VERTICAL:
            new_graph = new EntapGraphLineVertical(graphingData);
            break;

        case ENT_GRAPH_PIE_CHART:
            new_graph = new EntapGraphPieChart(graphingData);
            break;

        case ENT_GRAPH_BOX_PLOT_VERTICAL:
            new_graph = new EntapGraphBoxPlotVertical(graphingData);
            break;

        case ENT_GRAPH_BOX_PLOT_HORIZONTAL:
            new_graph = new EntapGraphBoxPlotHorizontal(graphingData);
            break;

        case ENT_GRAPH_BAR_STACKED:
            new_graph = new EntapGraphBarStacked(graphingData);
            break;

        default:
            ret = false;
            break;
    }

    if (new_graph != nullptr) {
        mGraphData.emplace(graphingData.text_file_path, new_graph);
    } else {
        ret = false;
    }

    return ret;
}

bool GraphingManager::add_datapoint(std::string &path, std::list<std::string> list) {
    bool ret = true;

    if (mGraphData.find(path) != mGraphData.end()) {
        mGraphData[path]->add_datapoint(list);
    } else {
        ret = false;
    }
    return ret;
}

void GraphingManager::graph_data(std::string &path) {
    EntapGraphBase *graph_base=nullptr;

    // If graphing enabled
    if (mGraphingEnabled) {
        // Yes, graphing is enabled
        if (mGraphData.find(path) != mGraphData.end()) {
            graph_base = mGraphData[path];
        }

        if (graph_base != nullptr) {
            graph(&graph_base->getMGraphingData());
            graph_base->close_graphing_file();
            SAFE_DELETE(graph_base);
            mGraphData.erase(path);
        }
    }
}

void GraphingManager::close_graphing_file(std::string &path) {

    EntapGraphBase *graph_base=nullptr;

    // If graphing enabled
    if (mGraphingEnabled) {
        // Yes, graphing is enabled
        if (mGraphData.find(path) != mGraphData.end()) {
            graph_base = mGraphData[path];
        }

        if (graph_base != nullptr) {
            graph_base->close_graphing_file();
            SAFE_DELETE(graph_base);
            mGraphData.erase(path);
        }
    }
}

//**********************************************************************
//**********************************************************************
//                   EntapGraphBase Nested Class
//**********************************************************************
//**********************************************************************

GraphingManager::EntapGraphBase::EntapGraphBase(GraphingManager::GraphingData &graphingData) {
    this->mGraphingData = graphingData;
    this->mpOutputTextStream = new std::ofstream(graphingData.text_file_path, std::ios::out | std::ios::trunc);
    *this->mpOutputTextStream << graphingData.x_axis_label << '\t' << graphingData.y_axis_label << std::endl;
}

GraphingManager::EntapGraphBase::~EntapGraphBase() {
    (void) this->close_graphing_file();
}


bool GraphingManager::EntapGraphBase::close_graphing_file() {
    bool ret = true;

    FS_dprint("GraphingManager - closing graphing file at: " + this->mGraphingData.text_file_path);

    if (mpOutputTextStream == nullptr) return ret; // File was already cleaned up
    try {
        if (this->mpOutputTextStream->is_open()) {
            this->mpOutputTextStream->close();
            SAFE_DELETE(this->mpOutputTextStream);
        }
    } catch (...) {
        SAFE_DELETE(this->mpOutputTextStream);
        ret = false;
    }
    return ret;
}

GraphingManager::GraphingData &GraphingManager::EntapGraphBase::getMGraphingData() {
    return this->mGraphingData;
}

//**********************************************************************
//**********************************************************************
//                      EntapGraphBarVertical Nested Class
//**********************************************************************
//**********************************************************************

GraphingManager::EntapGraphBarVertical::EntapGraphBarVertical(GraphingManager::GraphingData &graphingData) : EntapGraphBase(
        graphingData) {
}

bool GraphingManager::EntapGraphBarVertical::add_datapoint(std::list<std::string> &list) {
    if (list.size() < MINIMUM_ARGS) {
        return false;
    } else {
        *this->mpOutputTextStream << list.front() << '\t';
        list.pop_front();
        *this->mpOutputTextStream << list.front() << std::endl;
        return true;
    }
}

//**********************************************************************
//**********************************************************************
//                      EntapGraphBarHorizontal Nested Class
//**********************************************************************
//**********************************************************************

GraphingManager::EntapGraphBarHorizontal::EntapGraphBarHorizontal(GraphingManager::GraphingData &graphingData) : EntapGraphBase(
        graphingData) {
}

bool GraphingManager::EntapGraphBarHorizontal::add_datapoint(std::list<std::string> &list) {
    if (list.size() < MINIMUM_ARGS) {
        return false;
    } else {
        *this->mpOutputTextStream << list.front() << '\t';
        list.pop_front();
        *this->mpOutputTextStream << list.front() << std::endl;
        return true;
    }
}

//**********************************************************************
//**********************************************************************
//                 EntapGraphLineHorizontal Nested Class
//**********************************************************************
//**********************************************************************

GraphingManager::EntapGraphLineHorizontal::EntapGraphLineHorizontal(GraphingManager::GraphingData &graphingData)
        : EntapGraphBase(graphingData) {

}

bool GraphingManager::EntapGraphLineHorizontal::add_datapoint(std::list<std::string> &list) {

    return false;
}

//**********************************************************************
//**********************************************************************
//                 EntapGraphLineVertical Nested Class
//**********************************************************************
//**********************************************************************

GraphingManager::EntapGraphLineVertical::EntapGraphLineVertical(GraphingManager::GraphingData &graphingData)
        : EntapGraphBase(graphingData) {

}

bool GraphingManager::EntapGraphLineVertical::add_datapoint(std::list<std::string> &list) {
    return false;
}

//**********************************************************************
//**********************************************************************
//                 EntapGraphPieChart Nested Class
//**********************************************************************
//**********************************************************************

GraphingManager::EntapGraphPieChart::EntapGraphPieChart(GraphingManager::GraphingData &graphingData) : EntapGraphBase(
        graphingData) {

}

bool GraphingManager::EntapGraphPieChart::add_datapoint(std::list<std::string> &list) {
    bool ret = true;

    if (list.size() < MINIMUM_ARGS || list.size() > MAXIMUM_ARGS) {
        ret = false;
    } else {
        *mpOutputTextStream << list.front() << '\t';
        list.pop_front();
        *mpOutputTextStream << list.front() << std::endl;
    }
    return ret;
}

//**********************************************************************
//**********************************************************************
//                 EntapGraphBoxPlotVertical Nested Class
//**********************************************************************
//**********************************************************************

GraphingManager::EntapGraphBoxPlotVertical::EntapGraphBoxPlotVertical(GraphingManager::GraphingData &graphingData)
        : EntapGraphBase(graphingData) {

}

bool GraphingManager::EntapGraphBoxPlotVertical::add_datapoint(std::list<std::string> &list) {
    bool ret = true;

    if (list.size() < MINIMUM_ARGS || list.size() > MAXIMUM_ARGS) {
        ret = false;
    } else {
        *mpOutputTextStream << list.front() << '\t';
        list.pop_front();
        *mpOutputTextStream << list.front() << std::endl;
    }
    return ret;}

//**********************************************************************
//**********************************************************************
//                 EntapGraphBoxPlotHorizontal Nested Class
//**********************************************************************
//**********************************************************************

GraphingManager::EntapGraphBoxPlotHorizontal::EntapGraphBoxPlotHorizontal(GraphingManager::GraphingData &graphingData)
        : EntapGraphBase(graphingData) {

}

bool GraphingManager::EntapGraphBoxPlotHorizontal::add_datapoint(std::list<std::string> &list) {
    bool ret = true;

    if (list.size() < MINIMUM_ARGS || list.size() > MAXIMUM_ARGS) {
        ret = false;
    } else {
        *mpOutputTextStream << list.front() << '\t';
        list.pop_front();
        *mpOutputTextStream << list.front() << std::endl;
    }
    return ret;
}

//**********************************************************************
//**********************************************************************
//                 EntapGraphBarStacked Nested Class
//**********************************************************************
//**********************************************************************

GraphingManager::EntapGraphBarStacked::EntapGraphBarStacked(GraphingManager::GraphingData &graphingData)
        : EntapGraphBase(graphingData) {

}

bool GraphingManager::EntapGraphBarStacked::add_datapoint(std::list<std::string> &list) {
    bool ret = true;
    std::string out;

    if (list.size() < MINIMUM_ARGS || list.size() > MAXIMUM_ARGS) {
        ret = false;
    } else {
        for (std::string &str : list) {
            out += str + "\t";
        }

        out.pop_back(); // Remove trailing tab character
        *mpOutputTextStream << out << std::endl;
    }
    return ret;
}
