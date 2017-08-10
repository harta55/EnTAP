//
// Created by harta on 5/7/17.
//

#ifndef ENTAP_EXPRESSIONANALYSIS_H
#define ENTAP_EXPRESSIONANALYSIS_H
#include <iostream>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"
#include "GraphingManager.h"


class ExpressionAnalysis {
public:
    ExpressionAnalysis(std::string&, int, std::string&, std::string&,
                       boost::program_options::variables_map&, GraphingManager*);
    std::string execute(std::string, std::map<std::string, QuerySequence>& );

private:

    const std::string RSEM_OUT_DIR          = "expression/";
    const std::string RSEM_PROCESSED_DIR    = "processed/";
    const std::string RSEM_FIGURE_DIR       = "/figures";
    const std::string GRAPH_TXT_BOX_PLOT    = "comparison_box.txt";
    const std::string GRAPH_PNG_BOX_PLOT    = "comparison_box.png";
    const std::string GRAPH_TITLE_BOX_PLOT  = "Expression_Analysis";
    const std::string GRAPH_REJECTED_FLAG   = "Removed";
    const std::string GRAPH_KEPT_FLAG       = "Selected";
    const float REJECTED_ERROR_CUTOFF       = 75.0;

    const unsigned char GRAPH_EXPRESSION_FLAG = 2;
    const unsigned char GRAPH_BOX_FLAG        = 1;
    static constexpr int RSEM_COL_NUM = 7;

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

    std::string rsem(std::map<std::string, QuerySequence>&);
    std::string rsem_filter(std::string&,std::map<std::string, QuerySequence>&);
    bool is_file_empty(std::string);
    bool rsem_validate_file(std::string);
    bool rsem_conv_to_bam(std::string);
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
