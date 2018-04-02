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

#ifndef ENTAP_GRAPHINGMANAGER_H
#define ENTAP_GRAPHINGMANAGER_H
#include <iostream>
#include "EntapGlobals.h"

struct  GraphingData{
    std::string     text_file_path;
    std::string     graph_title;
    std::string     fig_out_path;
    uint8           software_flag;
    uint8           graph_type;
};


class GraphingManager {

public:
    GraphingManager(std::string);
    void graph(GraphingData&);
    bool is_graphing_enabled() const;

private:

    const std::string FLAG_GRAPH_TEXT   = "-i";
    const std::string FLAG_SOFT         = "-s";
    const std::string FLAG_GRAPH        = "-g";
    const std::string FLAG_TITLE        = "-t";
    const std::string FLAG_OUT_PATH     = "-p";

    std::string _graph_path;
    bool _graphing_enabled;


};


#endif //ENTAP_GRAPHING_H
