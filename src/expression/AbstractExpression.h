//
// Created by harta on 8/11/17.
//

#ifndef ENTAP_ABSTRACTEXPRESSION_H
#define ENTAP_ABSTRACTEXPRESSION_H

#include <string>
#include "../GraphingManager.h"
#include "../QuerySequence.h"

class AbstractExpression {
public:
    AbstractExpression(std::string &exe, std::string &out, std::string &in, std::string &proc,
                  std::string &fig, std::string &exp, std::string &align,GraphingManager *graphing){
        _exe_path = exe;
        _outpath = out;
        _inpath = in;
        _processed_path = proc;
        _figure_path = fig;
        _expression_outpath = exp;
        pGraphingManager = graphing;
        _alignpath = align;
    }

    virtual ~AbstractExpression() = default;
    virtual std::pair<bool, std::string> verify_files()=0;
    virtual void execute(std::map<std::string, QuerySequence>&) = 0;
    virtual std::string filter(std::map<std::string, QuerySequence>&) = 0;
    virtual void set_data(int, float, bool)=0;


protected:
    std::string _alignpath;
    std::string _exe_path;
    std::string _outpath;
    std::string _inpath;
    std::string _processed_path;
    std::string _figure_path;
    std::string _expression_outpath;
    GraphingManager *pGraphingManager;
};

#endif //ENTAP_ABSTRACTEXPRESSION_H
