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

#include "ModHorizontalGeneTransferDiamond.h"
#include "../ExceptionHandler.h"

std::vector<ENTAP_HEADERS> ModHorizontalGeneTransferDiamond::DEFAULT_HEADERS = {
        ENTAP_HEADER_HORIZONTALLY_TRANSFERRED_GENE
};

ModHorizontalGeneTransferDiamond::ModHorizontalGeneTransferDiamond(std::string &execution_stage_path, std::string &fasta_path,
                                                                   EntapDataPtrs &entap_data)
: AbstractHorizontalGeneTransfer(execution_stage_path, fasta_path, entap_data, "Horizontal Gene Transfer", DEFAULT_HEADERS){

    mSoftwareFlag = HGT_DIAMOND;
    mModDiamondDonor = std::unique_ptr<AbstractSimilaritySearch>(new ModDiamond(
            execution_stage_path, fasta_path, entap_data, mDonorDatabasePaths, "DIAMOND_Donor", HGT_DIAMOND_DONOR
    ));

    mModDiamondRecipient = std::unique_ptr<AbstractSimilaritySearch>(new ModDiamond(
            execution_stage_path, fasta_path, entap_data, mRecipientDatabasePaths, "DIAMOND_Recipient", HGT_DIAMOND_RECIPIENT
    ));
}

EntapModule::ModVerifyData ModHorizontalGeneTransferDiamond::verify_files() {
    EntapModule::ModVerifyData donor_databases_verified;
    EntapModule::ModVerifyData recipient_databases_verified;
    EntapModule::ModVerifyData return_databases_verified;

    donor_databases_verified = mModDiamondDonor->verify_files();
    recipient_databases_verified = mModDiamondRecipient->verify_files();

    return_databases_verified.files_exist = donor_databases_verified.files_exist && recipient_databases_verified.files_exist;

    return return_databases_verified;
}

/***
 * Horizontal Gene Transfer
 * This is a process of determining if a gene has been horizontally transferred.
 *  There are several stages to be performed:
 *
 *  1. Run Similarity Search with our transcriptome against all of our Donor and Recipient Databases
 *  2. Leverage alignment information to determine which genes are 'HGT candidates'. These 'may' be HGT genes, but
 *      we'll need to use more information to be sure. This is where the GFF file comes into play
 *  3. Using GFF input from user, determine upstream and downstream genes to our 'HGT candidates'
 *  4. If upstream/downstream genes match the taxonomy of our Candidate then this is a true HGT gene
 */
void ModHorizontalGeneTransferDiamond::execute() {

    // 1. Similarity Search against Donor and Recipient databases
    FS_dprint("HGT Running Similarity Serach against Donor and Recipient Databases");
    try {
        mModDiamondDonor->execute();
        mModDiamondRecipient->execute();
    } catch (const ExceptionHandler &e) {
        throw e;
    }
    FS_dprint("HGT run against donor and recipients complete, determining HGT candidates...");
    // 2.

}

void ModHorizontalGeneTransferDiamond::parse() {

}

bool ModHorizontalGeneTransferDiamond::is_executable(std::string &exe) {
    return false;
}

bool ModHorizontalGeneTransferDiamond::set_version() {
    return false;
}
