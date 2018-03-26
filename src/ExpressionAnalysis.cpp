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


//*********************** Includes *****************************
#include <boost/filesystem.hpp>
#include <csv.h>
#include <boost/regex.hpp>
#include <iomanip>
#include "ExpressionAnalysis.h"
#include "ExceptionHandler.h"
#include "EntapConfig.h"
#include "EntapGlobals.h"
#include "EntapExecute.h"
#include "common.h"
#include "expression/ModRSEM.h"
#include "FileSystem.h"
//**************************************************************


/**
 * ======================================================================
 * Function ExpressionAnalysis::ExpressionAnalysis(std::string &input,int t,
 *                                                 std::string &out,
 *                                                 boost::program_options::variables_map& user_flags,
 *                                                 GraphingManager *graph, QueryData *queryData)
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
 * @param t             - Thread count
 * @param out           - Main outfiles directory path
 * @param user_flags    - Boost parsed user flags
 * @param graph         - Ptr to graphing manager
 * @param queryData     - Ptr to query data
 *
 * @return              - None
 *
 * =====================================================================
 */
ExpressionAnalysis::ExpressionAnalysis(std::string &input,
                                       GraphingManager *graph,
                                       QueryData *queryData,
                                       FileSystem *filesystem,
                                       UserInput *userinput) {
    FS_dprint("Spawn object - ExpressionAnalysis");

    _pQueryData       = queryData;
    _pFileSystem      = filesystem;
    _pUserInput       = userinput;
    _pGraphingManager = graph;
    _inpath           = input;

    if (queryData == nullptr || filesystem == nullptr || userinput == nullptr) {
        throw ExceptionHandler("Unable to allocate memory to ExpressionAnalysis.",
            ERR_ENTAP_MEM_ALLOC);
    }

    _software_flag = ENTAP_EXECUTE::EXP_FLAG_RSEM;
    _threads       = _pUserInput->get_supported_threads();
    _exepath       = RSEM_EXE_DIR;
    _outpath       = _pFileSystem->get_root_path();
    _trim          = _pUserInput->has_input(UInput::INPUT_FLAG_TRIM);
    _overwrite     = _pUserInput->has_input(UInput::INPUT_FLAG_OVERWRITE);
    _issingle      = _pUserInput->has_input(UInput::INPUT_FLAG_SINGLE_END);
    if (_pUserInput->has_input(UInput::INPUT_FLAG_ALIGN)) { // Will be true
        _alignpath = _pUserInput->get_user_input<std::string>(UInput::INPUT_FLAG_ALIGN);
    }
    _fpkm          = _pUserInput->get_user_input<fp32>(UInput::INPUT_FLAG_FPKM);
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
    std::pair<bool, std::string>        verify_pair;
    std::unique_ptr<AbstractExpression> ptr;

    _inpath = input;
    if (_overwrite) _pFileSystem->delete_dir(_rsem_dir);
    _pFileSystem->create_dir(_rsem_dir);

    try {
        ptr = spawn_object();
        ptr->set_data(_threads, _fpkm, _issingle);  // Will remove later
        verify_pair = ptr->verify_files();
        if (!verify_pair.first) ptr->execute();
        output = ptr->filter();
    } catch (const ExceptionHandler &e) {throw e;}
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
        case ENTAP_EXECUTE::EXP_FLAG_RSEM:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _exepath,
                    _outpath,
                    _inpath,
                    _rsem_dir,
                    _alignpath,
                    _pGraphingManager,
                    _pQueryData,
                    _pFileSystem
            ));
        default:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _exepath,
                    _outpath,
                    _inpath,
                    _rsem_dir,
                    _alignpath,
                    _pGraphingManager,
                    _pQueryData,
                    _pFileSystem
            ));
    }
}


