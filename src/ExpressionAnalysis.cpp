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
ExpressionAnalysis::ExpressionAnalysis(std::string &input,int t, std::string &out,
                                       boost::program_options::variables_map& user_flags,
                                       GraphingManager *graph, QueryData *queryData) {
    FS_dprint("Spawn object - ExpressionAnalysis");
    _software_flag = 0;
    _query_data    = queryData;
    _inpath        = input;
    _threads       = t;
    _exepath       = RSEM_EXE_DIR;
    _outpath       = out;
    _trim          = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_TRIM);
    _overwrite     = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _issingle      = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_SINGLE_END);
    if (user_flags.count(ENTAP_CONFIG::INPUT_FLAG_ALIGN)) {
        _alignpath = user_flags[ENTAP_CONFIG::INPUT_FLAG_ALIGN].as<std::string>();
    }
    _fpkm          = user_flags[ENTAP_CONFIG::INPUT_FLAG_FPKM].as<fp32>();
    _rsem_dir      = PATHS(out, RSEM_OUT_DIR);
    _proc_dir      = PATHS(_rsem_dir, RSEM_PROCESSED_DIR);
    _figure_dir    = PATHS(_proc_dir,RSEM_FIGURE_DIR);
    _graphingManager = graph;
    SOFTWARE       = static_cast<ExpressionSoftware>(_software_flag);
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
    if (_overwrite) boostFS::remove_all(_rsem_dir);
    FS_create_dir(_rsem_dir);
    FS_create_dir(_figure_dir);
    FS_create_dir(_proc_dir);
    try {
        ptr = spawn_object();
        ptr->set_data(_threads, _fpkm, _issingle);
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
    switch (SOFTWARE) {
        case RSEM:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _exepath, _outpath, _inpath, _proc_dir, _figure_dir,
                    _rsem_dir, _alignpath, _graphingManager, _query_data
            ));
        default:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _exepath, _outpath, _inpath, _proc_dir, _figure_dir,
                    _rsem_dir, _alignpath, _graphingManager, _query_data
            ));
    }
}


