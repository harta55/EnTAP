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


#ifndef ENTAP_MODRSEM_H
#define ENTAP_MODRSEM_H


#include "AbstractExpression.h"
#include "../QueryData.h"

class ModRSEM : public AbstractExpression{

public:
    ModRSEM(std::string &exe, std::string &out, std::string &in, std::string &proc,
            std::string &fig, std::string &exp, std::string &align,GraphingManager *graphing,
            QueryData *query) :
    AbstractExpression(exe, out, in, proc, fig, exp, align,graphing, query){}

    virtual std::pair<bool, std::string> verify_files() override ;
    virtual void execute() override ;
    virtual std::string filter() override ;
    virtual void set_data(int, float, bool) override    ;

private:

    const std::string RSEM_PREP_REF_EXE     = "rsem-prepare-reference";
    const std::string RSEM_CALC_EXP_EXE     = "rsem-calculate-expression";
    const std::string RSEM_OUT_KEPT         = "_kept.fasta";
    const std::string RSEM_OUT_REMOVED      = "_removed.fasta";
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
    bool        _ispaired;

    bool rsem_validate_file(std::string);
    bool rsem_conv_to_bam(std::string);
    bool is_file_empty(std::string);
};


#endif //ENTAP_MODRSEM_H
