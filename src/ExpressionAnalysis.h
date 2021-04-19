/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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

#include "GraphingManager.h"
#include "QueryData.h"
#include "ExceptionHandler.h"
#include "EntapConfig.h"
#include "common.h"
#include "FileSystem.h"
#include "expression/AbstractExpression.h"
#include "expression/ModRSEM.h"

class AbstractExpression;

enum EXPRESSION_SOFTWARE {
    EXP_RSEM,
    EXP_COUNT
};


class ExpressionAnalysis {
    
public:
    ExpressionAnalysis(std::string&,EntapDataPtrs&);
    std::string execute(std::string);

private:

    const std::string RSEM_OUT_DIR          = "expression/";    // Name of root expression analysis directory

    std::string         mInFastaPath;       // FASTA path input from user to perform filtering on
                                            // *** MUST **** be original to maintain consistent headers
    std::string         mAlignPath;         // Absolute path to alignment file (BAM/SAM)
    std::string         mExpressionDir;     // Absolute path to root Expression Analysis directory
    bool                mOverwrite;         // TRUE if user would like to overwrite previous data
    uint16              mSoftwareFlag;      // FLag indicating Expression Analysis software module
    QueryData          *mpQueryData;        // Pointer to all Query data
    FileSystem         *mpFileSystem;       // Pointer to EnTAP filesystem
    UserInput          *mpUserInput;        // Pointer to User Input data
    EntapDataPtrs       mpEntapPtrs;        // All required data for execution

    std::unique_ptr<AbstractExpression> spawn_object();
};


#endif //ENTAP_EXPRESSIONANALYSIS_H
