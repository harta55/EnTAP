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

#ifndef ENTAP_HORIZONTALGENETRANSFER_H
#define ENTAP_HORIZONTALGENETRANSFER_H

#include "common.h"
#include "EntapGlobals.h"
#include "FileSystem.h"
#include "UserInput.h"
#include "horizontal_gene_transfer/AbstractHorizontalGeneTransfer.h"


class HorizontalGeneTransfer {
public:
    HorizontalGeneTransfer(std::string &transcriptome_path, EntapDataPtrs &entap_data);
    void execute();

private:
    std::string mTranscriptomePath;
    FileSystem    *mpFileSystem;
    UserInput     *mpUserInput;
    EntapDataPtrs       *mpEntapData;       // Pointer to all data required for execution
    std::string         mHGTDirectory;
    HORIZONTAL_GENE_TRANSFER_SOFTWARE mSoftwareFlag;

    const std::string HORIZONTAL_GENE_TRANSFER_DIR = "horizontal_gene_transfer";


    std::unique_ptr<AbstractHorizontalGeneTransfer> spawn_object();
};


#endif //ENTAP_HORIZONTALGENETRANSFER_H
