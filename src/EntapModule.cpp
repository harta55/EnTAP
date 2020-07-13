/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2020, Alexander Hart, Dr. Jill Wegrzyn
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
#include "QueryData.h"

std::vector<ENTAP_HEADERS> EntapModule::mModuleHeaders = {
        ENTAP_HEADER_QUERY
};

/**
 * ======================================================================
 * @class EntapModule
 *
 * Description          - Class for each program/module EnTAP supports
 *                      - Modules are responsible for Execution, Graphing,
 *                        calculation of statistics, verification of
 *                        required inputs to ensure execution is valid
 *
 * Notes                - None
 *
 * =====================================================================
 */
EntapModule::EntapModule(std::string &execution_stage_path, std::string &in_hits, EntapDataPtrs &entap_data,
                         std::string module_name, std::vector<ENTAP_HEADERS> &module_headers) {

    mOutpath = execution_stage_path;       // Should already be created
    mInputTranscriptome  = in_hits;

    mpGraphingManager = entap_data.mpGraphingManager;
    mpQueryData       = entap_data.mpQueryData;
    mpFileSystem      = entap_data.mpFileSystem;
    mpUserInput       = entap_data.mpUserInput;
    mpEntapDatabase   = entap_data.mpEntapDatabase;
    mModuleHeaders    = module_headers;

    // Initialize version numbers
    mVersionMajor = 0;
    mVersionMinor = 0;
    mVersionRev   = 0;

    mThreads         = mpUserInput->get_supported_threads();
    mBlastp          = mpUserInput->has_input(INPUT_FLAG_RUNPROTEIN);
    mOverwrite       = mpUserInput->has_input(INPUT_FLAG_OVERWRITE);
    mAlignmentFileTypes = mpUserInput->get_user_output_types();   // may be overridden at lower level
    mGoLevels        = mpUserInput->get_user_input<ent_input_multi_int_t >(INPUT_FLAG_GO_LEVELS);
    mEntapHeaders    = mpUserInput->get_user_input<std::vector<ENTAP_HEADERS>>(INPUT_FLAG_ENTAP_HEADERS);

    mModuleName      = module_name;
    mTranscriptomeShortname = mpFileSystem->get_filename(mInputTranscriptome, false);

    // INIT directories
    mModOutDir  = PATHS(mOutpath, module_name);
    mFigureDir  = PATHS(mModOutDir, FIGURE_DIR);
    mProcDir    = PATHS(mModOutDir, PROCESSED_OUT_DIR);
    mOverallResultsDir = PATHS(mModOutDir, OVERALL_RESULTS_DIR);    // generated at app level (some don't need this directory)

    // If overwriting data, remove entire execution stage directory
    if (mOverwrite) {
        mpFileSystem->delete_dir(mModOutDir);
    } else {
        mpFileSystem->delete_dir(mFigureDir);
        mpFileSystem->delete_dir(mProcDir);
    }
    // Recreate module + figure + processed directories
    mpFileSystem->create_dir(mModOutDir);
    mpFileSystem->create_dir(mFigureDir);
    mpFileSystem->create_dir(mProcDir);

    enable_headers();
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

void EntapModule::enable_headers() {
    for (auto &header : mModuleHeaders) {
        mpQueryData->header_set(header, true);
    }
}
