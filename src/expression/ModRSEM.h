/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
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


#ifndef ENTAP_MODRSEM_H
#define ENTAP_MODRSEM_H

//*********************** Includes *****************************
#include "AbstractExpression.h"
#include <csv.h>
#include "../ExceptionHandler.h"
#include "../GraphingManager.h"
#include "../FileSystem.h"
#include "../QuerySequence.h"
#include "../common.h"

//**************************************************************

class ModRSEM : public AbstractExpression {

public:
    ModRSEM(std::string &execution_stage_path, std::string &in_hits,
            EntapDataPtrs &entap_data, std::string &exe,
            std::string &align);

    ~ModRSEM() ;

    virtual std::pair<bool, std::string> verify_files() override ;
    virtual void execute() override ;
    virtual void parse() override;
    virtual void set_data(int, float, bool) override    ;

    virtual std::string get_final_fasta() override ;

private:

    const std::string RSEM_NAME             = "RSEM";

    const std::string RSEM_SAM_VALID        = "rsem-sam-validator";
    const std::string RSEM_PREP_REF_EXE     = "rsem-prepare-reference";
    const std::string RSEM_CALC_EXP_EXE     = "rsem-calculate-expression";
    const std::string RSEM_CONV_SAM         = "convert-sam-for-rsem";
    const std::string RSEM_OUT_KEPT         = "_kept.fasta";
    const std::string RSEM_OUT_REMOVED      = "_removed.fasta";
    const std::string RSEM_OUT_FILE         = ".genes.results";
    const std::string STD_REF_OUT           = "_rsem_reference";
    const std::string STD_EXP_OUT           = "_rsem_exp";
    const std::string STD_VALID_OUT         = "_rsem_validate";
    const std::string STD_CONVERT_SAM       = "_rsem_convert";
    const std::string GRAPH_TXT_BOX_PLOT    = "comparison_box.txt";
    const std::string GRAPH_PNG_BOX_PLOT    = "comparison_box.png";
    const std::string GRAPH_TITLE_BOX_PLOT  = "Expression_Analysis";
    const std::string GRAPH_REJECTED_FLAG   = "Removed";
    const std::string GRAPH_KEPT_FLAG       = "Selected";
    const float REJECTED_ERROR_CUTOFF       = 75.0;

    const unsigned char GRAPH_EXPRESSION_FLAG = 2;
    const unsigned char GRAPH_BOX_FLAG        = 1;
    static constexpr int RSEM_COL_NUM = 7;

    std::string _filename;
    std::string _rsem_out;
    std::string _exp_out;
    int         _threads;
    float       _fpkm;
    bool        _issingle;

    bool rsem_validate_file(std::string);
#if 0
    bool rsem_conv_to_bam(std::string); // Currently unused
#endif
    bool rsem_generate_reference(std::string&);
    bool rsem_expression_analysis(std::string&, std::string&);
};


#endif //ENTAP_MODRSEM_H
