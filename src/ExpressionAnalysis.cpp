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

#include <boost/filesystem.hpp>
#include <csv.h>
#include <boost/regex.hpp>
#include <iomanip>
#include "ExpressionAnalysis.h"
#include "ExceptionHandler.h"
#include "EntapConfig.h"
#include "EntapGlobals.h"
#include "EntapExecute.h"
#include "expression/ModRSEM.h"

namespace boostFS = boost::filesystem;

ExpressionAnalysis::ExpressionAnalysis(std::string &input,int t, std::string &out
    , boost::program_options::variables_map& user_flags, GraphingManager *graph) {
    print_debug("Spawn object - ExpressionAnalysis");
    _inpath = input;
    _threads = t;
    _exepath = RSEM_EXE_DIR;
    _outpath = out;
    _software_flag = 0;
    _overwrite = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _ispaired = (bool) user_flags.count("paired-end");
    if (user_flags.count(ENTAP_CONFIG::INPUT_FLAG_ALIGN)) {
        _alignpath = user_flags[ENTAP_CONFIG::INPUT_FLAG_ALIGN].as<std::string>();
    }
    _fpkm = user_flags[ENTAP_CONFIG::INPUT_FLAG_FPKM].as<float>();
    _rsem_dir = PATHS(out, RSEM_OUT_DIR);
    _proc_dir = PATHS(_rsem_dir, RSEM_PROCESSED_DIR);
    _figure_dir = PATHS(_proc_dir,RSEM_FIGURE_DIR);
    _graphingManager = graph;
    SOFTWARE = static_cast<ExpressionSoftware>(_software_flag);
}

std::string ExpressionAnalysis::execute(std::string input,
                                        std::map<std::string, QuerySequence>& MAP) {
    std::string output;
    std::pair<bool, std::string> verify_pair;
    std::unique_ptr<AbstractExpression> ptr;

    _inpath = input;
    if (_overwrite) boostFS::remove_all(_rsem_dir);
    boostFS::create_directories(_rsem_dir);
    boostFS::create_directories(_figure_dir);
    boostFS::create_directories(_proc_dir);
    try {
        ptr = spawn_object();
        ptr->set_data(_threads, _fpkm, _ispaired);
        verify_pair = ptr->verify_files();
        if (!verify_pair.first) ptr->execute(MAP);
        output = ptr->filter(MAP);
        ptr.release();
    } catch (const ExceptionHandler &e) {throw e;}
    return output;
}

std::unique_ptr<AbstractExpression> ExpressionAnalysis::spawn_object() {
    switch (SOFTWARE) {
        case RSEM:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _exepath, _outpath, _inpath, _proc_dir, _figure_dir,
                    _rsem_dir, _alignpath, _graphingManager
            ));
        default:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _exepath, _outpath, _inpath, _proc_dir, _figure_dir,
                    _rsem_dir, _alignpath, _graphingManager
            ));
    }
}


