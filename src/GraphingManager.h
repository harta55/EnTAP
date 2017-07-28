//
// Created by harta on 7/21/17.
//

#ifndef ENTAP_GRAPHINGMANAGER_H
#define ENTAP_GRAPHINGMANAGER_H
#include <iostream>
#include "EntapConsts.h"

class GraphingManager {


public:
    GraphingManager(std::string);
    void graph(GraphingStruct&);


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
