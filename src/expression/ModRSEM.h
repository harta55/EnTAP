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

/**
 * ======================================================================
 * @class ModRSEM
 *
 * Description          - This EnTAP module supports execution, parsing, and
 *                        statistical analysis of the RSEM software
 *                        through terminal commands
 *                      - RSEM performs Expression Analysis through processing
 *                        reads from an alignment file input by the user
 *                      - Parsed data is added to QueryData class
 *                      - Transcriptome will be filtered based on FPKM value
 *                      - Inherits from AbstractExpression and EntapModule classes
 *
 * Citation             - B. Li and C. N. Dewey, “RSEM: accurate transcript
 *                        quantification from RNA-Seq data with or without
 *                        a reference genome,” (in eng), BMC Bioinformatics,
 *                        vol. 12, p. 323, Aug 2011.
 *
 * ======================================================================
 */
class ModRSEM : public AbstractExpression {

public:
    ModRSEM(std::string &execution_stage_path, std::string &in_hits,
            EntapDataPtrs &entap_data, std::string &align);

    ~ModRSEM();

    virtual ModVerifyData verify_files() override ;
    virtual void execute() override ;
    virtual void parse() override;
    virtual std::string get_final_fasta() override ;
    virtual bool set_version() override;

private:

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

    static constexpr int RSEM_COL_NUM = 7;
    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;

    std::string mFilename;
    std::string mRsemOut;
    std::string mExpressionOut;
    std::string mCalcExpressionExe;     // User input to RSEM calc expression exe
    std::string mSamValidExe;           // User input to RSEM sam validate exe
    std::string mPrepReferenceExe;      // User input to RSEM prep reference exe
    std::string mConvertSamExe;         // User input to RSEM convert SAM for RSEM exe
    bool        mIsSingle;

    bool rsem_validate_file(std::string);
#if 0
    bool rsem_conv_to_bam(std::string); // Currently unused
#endif
    bool rsem_generate_reference(std::string&);
    bool rsem_expression_analysis(std::string&, std::string&);
};


#endif //ENTAP_MODRSEM_H
