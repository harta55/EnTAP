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


#ifndef ENTAP_ABSTRACTEXPRESSION_H
#define ENTAP_ABSTRACTEXPRESSION_H

//*********************** Includes *****************************
#include "../common.h"
#include "../EntapModule.h"
#include "../GraphingManager.h"

//**************************************************************


/**
 * ======================================================================
 * Class AbstractExpression
 *
 * Description          - Abstract class/functions for Expression Analysis
 *                        portion of pipeline
 *                      - New modules will inherit from this class and
 *                        implement it's functions to fit into pipeline
 *
 * Notes                - None
 *
 *
 * @return              - None
 * ======================================================================
 */
class AbstractExpression : public EntapModule {

public:

/**
 * ======================================================================
 * Function AbstractExpression(std::string &exe, std::string &out,
 *                             std::string &in, std::string &proc,
                               std::string &fig, std::string &exp,
                               std::string &align,GraphingManager *graphing,
                               QueryData *query, bool trim)
 *
 * Description          - Constructor for Abstract expression class
 *                      - Initializes protected member variables for
 *                        expression modules
 *
 * Notes                - Constructor
 *
 * @param exe           - Path to execution directory (EnTAP, unused)
 * @param out           - Path to main outfiles directory (unused)
 * @param in            - Path to original transcriptome from user
 * @param proc          - Path to processed directory (within rsem for now)
 * @param fig           - Path to figure directory (within rsem for now)
 * @param exp           - Path to expression directory
 * @param align         - Path to SAM/BAM alignment file
 * @param graphing      - Ptr to graphing manager
 * @param query         - Pointer to query data
 * @param trim          - Bool of yes (sequences to be trimmed) or no
 *
 * @return              - None
 * ======================================================================
 */
    AbstractExpression(std::string &execution_stage_path, std::string &in_hits,
                       EntapDataPtrs &entap_data, std::string module_name,
                       std::vector<ENTAP_HEADERS> &module_headers);

    virtual ~AbstractExpression() = default;
    virtual ModVerifyData verify_files()=0;
    virtual void execute() = 0;
    virtual void parse() = 0;
    virtual std::string get_final_fasta()=0;
    virtual bool set_version() = 0;

    virtual void set_success_flags() override ;

protected:
    std::string     mAlignPath;     // Absolute path to alignment file (BAM/SAM)
    std::string     mFinalFasta;    // Absolute path to final filtered FASTA produced from this stage
    fp64            mFPKM;          // FPKM threshold user would like to filter by
    bool            mIsSingle;       // TRUE if use has non-paired data
};

#endif //ENTAP_ABSTRACTEXPRESSION_H
