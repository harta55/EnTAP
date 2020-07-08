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

#include "ModBUSCO.h"

std::vector<ENTAP_HEADERS> ModBUSCO::DEFAULT_HEADERS = {
    ENTAP_HEADER_ONT_BUSCO_ID,
    ENTAP_HEADER_ONT_BUSCO_STATUS,
    ENTAP_HEADER_ONT_BUSCO_LENGTH,
    ENTAP_HEADER_ONT_BUSCO_SCORE
};

ModBUSCO::ModBUSCO(std::string &in_hits, std::string &ont_out, EntapDataPtrs &entap_data) :
        AbstractOntology(in_hits, ont_out, entap_data, "BUSCO", DEFAULT_HEADERS) {
    FS_dprint("Spawn Object - ModBUSCO");

    // Pull inputs from user
    mExePath       = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_BUSCO_EXE);
    mDatabasePaths = mpUserInput->get_user_input<ent_input_multi_str_t>(INPUT_FLAG_BUSCO_DATABASE);
    mExePath       = mpUserInput->get_user_input<ent_input_fp_t >(INPUT_FLAG_BUSCO_EVAL);
    mBlastp ? mRunType = BUSCO_RUN_TYPE_PROT : mRunType= BUSCO_RUN_TYPE_TRAN;
}

ModBUSCO::~ModBUSCO() {
    FS_dprint("Killing object - ModBUSCO");
}

bool ModBUSCO::is_executable(std::string &exe) {
    std::string test_command;
    TerminalData terminalData;

    test_command = exe + " --help";

    terminalData.command = test_command;
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;

    return TC_execute_cmd(terminalData) == 0;
}

void ModBUSCO::execute() {
    std::string  database_out_path;     // Base path that BUSCO will output files to
                                        // WARNING!!! BUSCO prepends "run_" to this
    std::string  database_out_filename;
    std::string  std_out;
    std::string  cmd;
    TerminalData terminalData;

    // first we want to get the version

    // Execute BUSCO analysis against all of the databases input by the user
    for (std::string &database : mDatabasePaths) {

        if (mpFileSystem->file_exists(database)) {
            FS_dprint("Running BUSCO against database at: " + database);

            switch (mBuscoVersion) {

                case BUSCO_VERSION_3:
                    database_out_filename = mpFileSystem->get_filename(database, false) + "_" + mRunType + "_results";
                    database_out_path = PATHS(mModOutDir, database_out_filename);
                    std_out = PATHS(mModOutDir, database_out_filename) + "_" + FileSystem::EXT_STD;

                    cmd =
                            mExePath + " " +
                            BUSCO_INPUT_IN       + " " + mInputTranscriptome + " " + /* Input Transcriptome*/
                            BUSCO_INPUT_OUTPUT   + " " + database_out_path   + " " + /* Output base path */
                            BUSCO_INPUT_RUN_TYPE + " " + mRunType            + " " + /* Run type (tran or prot) */
                            BUSCO_INPUT_DATABASE + " " + database            + " " + /* Database */
                            BUSCO_INPUT_CPU      + " " + std::to_string(mThreads);

                    terminalData = {};
                    terminalData.base_std_path = std_out;
                    terminalData.command = cmd;
                    terminalData.print_files = true;
                    terminalData.suppress_std_err = false;

                    if (TC_execute_cmd(terminalData) != 0) {
                        // Error in run
                    }
                    break;

                case BUSCO_VERSION_4:
                    break;

                default:
                    break;
            }

        } else {
            // NO, Database does not exist!! Default exit execution
            throw ExceptionHandler("ERROR BUSCO database does not exist at: " + database,
                ERR_ENTAP_RUN_BUSCO);
        }
    }
}

void ModBUSCO::parse() {

}

EntapModule::ModVerifyData ModBUSCO::verify_files() {
    return EntapModule::ModVerifyData();
}

std::string ModBUSCO::get_output_dir(std::string &output_path) {
    std::string ret;

    return ret;
}

void ModBUSCO::get_version() {
    return;
}
