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
    std::string _inpath,_alignpath, _exepath, _outpath;
    bool _ispaired, _overwrite;
    int _threads;
    float _fpkm;
    short _software_flag;
    std::string rsem(std::map<std::string, QuerySequence>&);
    std::string rsem_filter(std::string&,std::map<std::string, QuerySequence>&);
    bool is_file_empty(std::string);
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
