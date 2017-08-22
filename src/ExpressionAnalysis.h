//
// Created by harta on 5/7/17.
//

#ifndef ENTAP_EXPRESSIONANALYSIS_H
#define ENTAP_EXPRESSIONANALYSIS_H
#include <iostream>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"
#include "GraphingManager.h"
#include "expression/AbstractExpression.h"


class ExpressionAnalysis {
public:
    ExpressionAnalysis(std::string&, int, std::string&,boost::program_options::variables_map&,
                       GraphingManager*);
    std::string execute(std::string, std::map<std::string, QuerySequence>& );

private:

    enum ExpressionSoftware {
        RSEM
    };


    const std::string RSEM_OUT_DIR          = "expression/";
    const std::string RSEM_PROCESSED_DIR    = "processed/";
    const std::string RSEM_FIGURE_DIR       = "/figures";

    std::string         _inpath;
    std::string         _alignpath;
    std::string         _exepath;
    std::string         _outpath;
    std::string         _rsem_dir;
    std::string         _proc_dir;
    std::string         _figure_dir;
    bool                _ispaired;
    bool                _overwrite;
    short               _software_flag;
    int                 _threads;
    float               _fpkm;
    GraphingManager  *_graphingManager;
    ExpressionSoftware SOFTWARE;

    std::unique_ptr<AbstractExpression> spawn_object();
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
