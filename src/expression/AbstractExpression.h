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
