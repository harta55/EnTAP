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

    std::string         _inpath;
    std::string         _alignpath;
    std::string         _exepath;
    std::string         _outpath;
    bool                _ispaired;
    bool                _overwrite;
    short               _software_flag;
    int                 _threads;
    float               _fpkm;

    std::string rsem(std::map<std::string, QuerySequence>&);
    std::string rsem_filter(std::string&,std::map<std::string, QuerySequence>&);
    bool is_file_empty(std::string);
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
