/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/


#ifndef ENTAP_ABSTRACTFRAME_H
#define ENTAP_ABSTRACTFRAME_H

#include "../QuerySequence.h"

class AbstractFrame {
public:
    AbstractFrame(std::string &exe, std::string &out, std::string &in, std::string &proc,
                std::string &fig, std::string &frame, GraphingManager *graphing){
        _exe_path = exe;
        _outpath = out;
        _inpath = in;
        _processed_path = proc;
        _figure_path = fig;
        _frame_outpath = frame;
        pGraphingManager = graphing;

    }
    virtual ~AbstractFrame() = default;
    virtual std::pair<bool, std::string> verify_files()=0;
    virtual std::string execute(std::map<std::string, QuerySequence>&) = 0;
    virtual void parse(std::map<std::string, QuerySequence>&) = 0;


protected:
    std::string _exe_path;
    std::string _outpath;
    std::string _inpath;
    std::string _processed_path;
    std::string _figure_path;
    std::string _frame_outpath;
    GraphingManager *pGraphingManager;
};


#endif //ENTAP_ABSTRACTFRAME_H
