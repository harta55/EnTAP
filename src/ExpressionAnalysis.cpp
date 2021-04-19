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


//*********************** Includes *****************************
#include "ExpressionAnalysis.h"

//**************************************************************


/**
 * ======================================================================
 * Function ExpressionAnalysis::ExpressionAnalysis(std::string &input, EntapDataPtrs entap_data)
 *
 * Description          - Main expression analysis constructor
 *                      - Handles execution of specific packages as well
 *                        as managing overall results (not currently)
 *                      - Initializes member variables that will be
 *                        pushed to specific modules
 *
 * Notes                - Constructor
 *
 * @param input         - Path to original input fasta file (may change)
 * @param entap_data    - Entap data pointers
 *
 * @return              - None
 *
 * =====================================================================
 */
ExpressionAnalysis::ExpressionAnalysis(std::string &input,EntapDataPtrs& entap_data) {
    FS_dprint("Spawn Object - ExpressionAnalysis");

    mpQueryData       = entap_data.mpQueryData;
    mpFileSystem      = entap_data.mpFileSystem;
    mpUserInput       = entap_data.mpUserInput;
    mInFastaPath           = input;

    mpEntapPtrs = entap_data;

    mSoftwareFlag = EXP_RSEM;
    mOverwrite     = mpUserInput->has_input(INPUT_FLAG_OVERWRITE);
    if (mpUserInput->has_input(INPUT_FLAG_ALIGN)) { // Will be true
        mAlignPath = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_ALIGN);
    }
    mExpressionDir      = PATHS(mpFileSystem->get_root_path(), RSEM_OUT_DIR);
}


/**
 * ======================================================================
 * Function std::string ExpressionAnalysis::execute(std::string input)
 *
 * Description          - Entry into expression analysis
 *                      - Spawns whichever module user selects (currently
 *                        only RSEM)
 *                      - Creates/removes output directories as instructed
 *
 * Notes                - Entry
 *
 * @param input         - Path to original input fasta file (may change)
 *
 * @return              - Path to filtered fasta file
 *
 * =====================================================================
 */
std::string ExpressionAnalysis::execute(std::string input) {
    std::string                         output;
    EntapModule::ModVerifyData          verify_data;
    std::unique_ptr<AbstractExpression> ptr;

    mInFastaPath = input;
    if (mOverwrite) mpFileSystem->delete_dir(mExpressionDir);
    mpFileSystem->create_dir(mExpressionDir);

    try {
        ptr = spawn_object();
        verify_data = ptr->verify_files();
        if (!verify_data.files_exist) ptr->execute();
        ptr->parse();
        output = ptr->get_final_fasta();

        // If successful, set flags
        ptr->set_success_flags();
        ptr.reset();

    } catch (const ExceptionHandler &e) {
        ptr.reset();
        throw e;
    }
    return output;
}


/**
 * ======================================================================
 * Function std::unique_ptr<AbstractExpression> ExpressionAnalysis::spawn_object()
 *
 * Description          - Spawns module to be used within pipeline (currently
 *                        only RSEM)
 *
 * Notes                - None
 *
 *
 * @return              - Unique ptr (c++11) for module
 *
 * =====================================================================
 */
std::unique_ptr<AbstractExpression> ExpressionAnalysis::spawn_object() {
    switch (mSoftwareFlag) {
        case EXP_RSEM:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    mExpressionDir,
                    mInFastaPath,
                    mpEntapPtrs,
                    mAlignPath
            ));
        default:
            FS_dprint("WARNING unexpected software type in ExpressionAnalysis");
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    mExpressionDir,
                    mInFastaPath,
                    mpEntapPtrs,
                    mAlignPath
            ));
    }
}


