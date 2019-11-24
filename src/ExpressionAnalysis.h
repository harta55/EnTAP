/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
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

#include "QuerySequence.h"
#include "GraphingManager.h"
#include "QueryData.h"
#include "ExceptionHandler.h"
#include "EntapConfig.h"
#include "EntapGlobals.h"
#include "common.h"
#include "FileSystem.h"
#include "expression/AbstractExpression.h"
#include "expression/ModRSEM.h"

class AbstractExpression;

class ExpressionAnalysis {
public:
    ExpressionAnalysis(std::string&,EntapDataPtrs&);
    std::string execute(std::string);

private:

    const std::string RSEM_OUT_DIR          = "expression/";

    std::string         _inpath;
    std::string         _alignpath;
    std::string         _exepath;
    std::string         _outpath;
    std::string         _rsem_dir;
    bool                _issingle;
    bool                _overwrite;
    short               _software_flag;
    int                 _threads;
    float               _fpkm;
    GraphingManager  *_pGraphingManager;
    QueryData        *_pQueryData;
    FileSystem       *_pFileSystem;
    UserInput        *_pUserInput;
    EntapDataPtrs    _entap_data;

    std::unique_ptr<AbstractExpression> spawn_object();
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
