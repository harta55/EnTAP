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

#include "HorizontalGeneTransfer.h"
#include "EntapModule.h"
#include "ExceptionHandler.h"
#include "horizontal_gene_transfer/ModHorizontalGeneTransferDiamond.h"

HorizontalGeneTransfer::HorizontalGeneTransfer(std::string &transcriptome_path, EntapDataPtrs &entap_data) {
    FS_dprint("Spawn Object - HorizontalGeneTransfer");
    mTranscriptomePath = transcriptome_path;
    mpFileSystem = entap_data.mpFileSystem;
    mpUserInput = entap_data.mpUserInput;
    mpEntapData = &entap_data;

    // Set HGT paths/directories
    mHGTDirectory  = PATHS(mpFileSystem->get_root_path(), HORIZONTAL_GENE_TRANSFER_DIR);

    if (mpUserInput->has_input(INPUT_FLAG_OVERWRITE)) {
        mpFileSystem->delete_dir(mHGTDirectory);
    }
    mpFileSystem->create_dir(mHGTDirectory);
    mSoftwareFlag = HGT_DIAMOND;
}

void HorizontalGeneTransfer::execute() {

    EntapModule::ModVerifyData verifyData;
    std::unique_ptr<AbstractHorizontalGeneTransfer> ptr;
    TC_print(TC_PRINT_COUT, "Beginning Horizontal Gene Transfer...");

    try {
        ptr = spawn_object();
        verifyData = ptr->verify_files();
        if (!verifyData.files_exist) {
            ptr->execute();
        }
        ptr->parse();
        ptr->set_success_flags();

        ptr.reset();
    } catch (const ExceptionHandler &e) {
        ptr.reset();
        throw e;
    }
    TC_print(TC_PRINT_COUT, "Horizontal Gene Transfer complete");

}

std::unique_ptr<AbstractHorizontalGeneTransfer> HorizontalGeneTransfer::spawn_object() {
    switch (mSoftwareFlag) {
        case HGT_DIAMOND:
        default:
            return std::unique_ptr<AbstractHorizontalGeneTransfer>(new ModHorizontalGeneTransferDiamond(
                    mHGTDirectory, mTranscriptomePath, *mpEntapData
            ));
    }
}
