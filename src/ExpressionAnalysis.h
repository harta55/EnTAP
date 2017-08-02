//
// Created by harta on 5/7/17.
//

#ifndef ENTAP_EXPRESSIONANALYSIS_H
#define ENTAP_EXPRESSIONANALYSIS_H
#include <iostream>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"


class ExpressionAnalysis {
public:
    ExpressionAnalysis(std::string&, int, std::string&, std::string&,
                       boost::program_options::variables_map&);
    std::string execute(std::string, std::map<std::string, QuerySequence>& );

private:

    const std::string RSEM_OUT_DIR          = "expression/";
    const std::string RSEM_PROCESSED_DIR    = "processed/";
    static constexpr int RSEM_COL_NUM = 7;

    std::string         _inpath;
    std::string         _alignpath;
    std::string         _exepath;
    std::string         _outpath;
    std::string         _rsem_dir;
    std::string         _proc_dir;
    bool                _ispaired;
    bool                _overwrite;
    short               _software_flag;
    int                 _threads;
    float               _fpkm;

    std::string rsem(std::map<std::string, QuerySequence>&);
    std::string rsem_filter(std::string&,std::map<std::string, QuerySequence>&);
    bool is_file_empty(std::string);
    bool rsem_validate_file(std::string);
    bool rsem_conv_to_bam(std::string);
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
