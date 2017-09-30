/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
 *
 * This file is part of EnTAP.
 *
 * EnTAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EnTAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ENTAP_EXPRESSIONANALYSIS_H
#define ENTAP_EXPRESSIONANALYSIS_H
#include <iostream>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"
#include "GraphingManager.h"
#include "expression/AbstractExpression.h"
#include "QueryData.h"


class ExpressionAnalysis {
public:
    ExpressionAnalysis(std::string&, int, std::string&,boost::program_options::variables_map&,
                       GraphingManager*, QueryData*);
    std::string execute(std::string);

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
    QueryData        *_query_data;
    ExpressionSoftware SOFTWARE;

    std::unique_ptr<AbstractExpression> spawn_object();
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
