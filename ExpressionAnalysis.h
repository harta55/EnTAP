//
// Created by harta on 5/7/17.
//

#ifndef ENTAP_EXPRESSIONANALYSIS_H
#define ENTAP_EXPRESSIONANALYSIS_H
#include <iostream>


class ExpressionAnalysis {
public:
    ExpressionAnalysis(std::string&, int, std::string&, std::string&, bool);
    std::string execute(short, bool, std::string);

private:
    std::string _inpath,_alignpath, _exepath, _outpath;
    bool _ispaired, _overwrite;
    int _threads;
    std::string rsem();
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
