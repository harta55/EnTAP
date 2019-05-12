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

#include "EntapModule.h"

EntapModule::EntapModule(std::string &execution_stage_path, std::string &in_hits, EntapDataPtrs &entap_data,
                         std::string module_name, std::string &exe_path) {

    _outpath  = execution_stage_path;       // Should already be created
    _in_hits  = in_hits;
    _exe_path = exe_path;

    _pGraphingManager = entap_data._pGraphingManager;
    _pQUERY_DATA      = entap_data._pQueryData;
    _pFileSystem      = entap_data._pFileSystem;
    _pUserInput       = entap_data._pUserInput;
    _pEntapDatabase   = entap_data._pEntapDatbase;

    _threads         = _pUserInput->get_supported_threads();
    _blastp          = _pUserInput->has_input(_pUserInput->INPUT_FLAG_RUNPROTEIN);
    _overwrite       = _pUserInput->has_input(_pUserInput->INPUT_FLAG_OVERWRITE);
    _alignment_file_types = _pUserInput->get_user_output_types();   // may be overridden at lower level

    _transcript_shortname = _pFileSystem->get_filename(_in_hits, false);

    // INIT directories
    _mod_out_dir = PATHS(_outpath, module_name);
    _figure_dir  = PATHS(_mod_out_dir, FIGURE_DIR);
    _proc_dir    = PATHS(_mod_out_dir, PROCESSED_OUT_DIR);
    _overall_results_dir = PATHS(_mod_out_dir, OVERALL_RESULTS_DIR);    // generated at app level

    // If overwriting data, remove entire execution stage directory
    if (_overwrite) {
        _pFileSystem->delete_dir(_mod_out_dir);
    } else {
        _pFileSystem->delete_dir(_figure_dir);
        _pFileSystem->delete_dir(_proc_dir);
    }
    // Recreate module + figure + processed directories
    _pFileSystem->create_dir(_mod_out_dir);
    _pFileSystem->create_dir(_figure_dir);
    _pFileSystem->create_dir(_proc_dir);
}


go_format_t EntapModule::EM_parse_go_list(std::string list, EntapDatabase* database,char delim) {

    go_format_t output;
    std::string temp;
    std::vector<std::vector<std::string>>results;

    if (list.empty()) return output;
    std::istringstream ss(list);
    while (std::getline(ss,temp,delim)) {
        GoEntry term_info =
                database->get_go_entry(temp);
        output[term_info.category].push_back(temp + "-" + term_info.term +
                                             "(L=" + term_info.level + ")");
    }
    return output;
}