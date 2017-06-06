//
// Created by harta on 5/7/17.
//

#ifndef ENTAP_EXPRESSIONANALYSIS_H
#define ENTAP_EXPRESSIONANALYSIS_H
#include <iostream>
#include <boost/program_options/variables_map.hpp>


class ExpressionAnalysis {
public:
    ExpressionAnalysis(std::string&, int, std::string&, std::string&,
                       boost::program_options::variables_map&);
    std::string execute();

private:
    std::string _inpath,_alignpath, _exepath, _outpath;
    bool _ispaired, _overwrite;
    int _threads;
    short _software_flag;
    std::string rsem();
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
