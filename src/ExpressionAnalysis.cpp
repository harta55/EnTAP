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

    _pQueryData       = entap_data._pQueryData;
    _pFileSystem      = entap_data._pFileSystem;
    _pUserInput       = entap_data._pUserInput;
    _pGraphingManager = entap_data._pGraphingManager;
    _inpath           = input;

    _entap_data = entap_data;

    _software_flag = EXP_RSEM;
    _threads       = _pUserInput->get_supported_threads();
    _exepath       = RSEM_EXE_DIR;
    _outpath       = _pFileSystem->get_root_path();
    _trim          = _pUserInput->has_input(_pUserInput->INPUT_FLAG_TRIM);
    _overwrite     = _pUserInput->has_input(_pUserInput->INPUT_FLAG_OVERWRITE);
    _issingle      = _pUserInput->has_input(_pUserInput->INPUT_FLAG_SINGLE_END);
    if (_pUserInput->has_input(_pUserInput->INPUT_FLAG_ALIGN)) { // Will be true
        _alignpath = _pUserInput->get_user_input<std::string>(_pUserInput->INPUT_FLAG_ALIGN);
    }
    _fpkm          = _pUserInput->get_user_input<fp32>(_pUserInput->INPUT_FLAG_FPKM);
    _rsem_dir      = PATHS(_outpath, RSEM_OUT_DIR);
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

    _inpath = input;
    if (_overwrite) _pFileSystem->delete_dir(_rsem_dir);
    _pFileSystem->create_dir(_rsem_dir);

    try {
        ptr = spawn_object();
        ptr->set_data(_threads, _fpkm, _issingle);  // Will remove later
        verify_data = ptr->verify_files();
        if (!verify_data.files_exist) ptr->execute();
        ptr->parse();
        output = ptr->get_final_fasta();
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
    switch (_software_flag) {
        case EXP_RSEM:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _rsem_dir,
                    _inpath,
                    _entap_data,
                    _exepath,
                    _alignpath
            ));
        default:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _rsem_dir,
                    _inpath,
                    _entap_data,
                    _exepath,
                    _alignpath
            ));
    }
}


