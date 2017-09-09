/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/


#ifndef ENTAP_FRAMESELECTION_H
#define ENTAP_FRAMESELECTION_H

//*********************** Includes *****************************
#include <iostream>
#include <map>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"
#include "GraphingManager.h"
#include "frame_selection/AbstractFrame.h"
//**************************************************************


class FrameSelection {


public:
    std::string execute(std::string,std::map<std::string,QuerySequence>&);
    FrameSelection(std::string&, std::string&, boost::program_options::variables_map &, GraphingManager*);


private:

    enum FrameSoftware {
        GENEMARKST
    };

    const std::string FRAME_SELECTION_OUT_DIR       = "frame_selection/";
    const std::string PROCESSED_DIR                 = "processed/";
    const std::string FIGURE_DIR                    = "figures/";

    std::string      _frame_outpath;
    std::string      _processed_path;
    std::string      _figure_path;
    std::string      _exe_path;
    std::string      _inpath;
    std::string      _outpath;
    bool             _overwrite;
    short            _software_flag;
    GraphingManager  *_graphingManager;
    FrameSoftware    SOFTWARE;

    std::unique_ptr<AbstractFrame> spawn_object();
};


#endif //ENTAP_FRAMESELECTION_H
