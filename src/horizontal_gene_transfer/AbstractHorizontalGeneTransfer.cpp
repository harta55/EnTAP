/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2023, Alexander Hart, Dr. Jill Wegrzyn
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

#include "AbstractHorizontalGeneTransfer.h"
#include "../EntapGlobals.h"
#include "../common.h"

AbstractHorizontalGeneTransfer::AbstractHorizontalGeneTransfer(std::string &execution_stage_path, std::string &in_hits,
                                                   EntapDataPtrs &entap_data, std::string mod_name, std::vector<ENTAP_HEADERS> &module_headers)
        :EntapModule(execution_stage_path, in_hits, entap_data, mod_name, module_headers) {

    mExecutionState = HORIZONTAL_GENE_TRANSFER;

    // Get relevant user info for Horizontal Gene Transfer
    mDonorDatabasePaths = mpUserInput->get_user_input<ent_input_multi_str_t>(INPUT_FLAG_HGT_DONOR_DATABASES);
    mRecipientDatabasePaths = mpUserInput->get_user_input<ent_input_multi_str_t>(INPUT_FLAG_HGT_RECIPIENT_DATABASES);
    mGFFPath = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_HGT_GFF);
    mQCoverage        = mpUserInput->get_user_input<ent_input_fp_t >(INPUT_FLAG_QCOVERAGE);
    mTCoverage        = mpUserInput->get_user_input<ent_input_fp_t >(INPUT_FLAG_TCOVERAGE);
    mEVal            = mpUserInput->get_user_input<ent_input_fp_t>(INPUT_FLAG_E_VALUE);

    // set blast string to use for file naming
    mBlastp ? mBlastType = BLASTP_STR : mBlastType = BLASTX_STR;

//    // create overall results dir
    mpFileSystem->delete_dir(mOverallResultsDir);
    mpFileSystem->create_dir(mOverallResultsDir);
}

void AbstractHorizontalGeneTransfer::set_success_flags() {
    mpQueryData->set_is_success_hgt(true);
}

std::string AbstractHorizontalGeneTransfer::get_database_shortname(std::string &full_path) {
    return mpFileSystem->get_filename(full_path, false);
}

std::string AbstractHorizontalGeneTransfer::get_database_output_path(std::string &database_name) {
    return PATHS(mModOutDir,mBlastType + "_" + mTranscriptomeShortname + "_" + get_database_shortname(database_name) + FileSystem::EXT_OUT);
}
