/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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
#include "FrameSelection.h"
#include "ExceptionHandler.h"
#include "frame_selection/ModGeneMarkST.h"
#include "FileSystem.h"
#include "frame_selection/ModTransdecoder.h"
//**************************************************************


/**
 * ======================================================================
 * Function FrameSelection(std::string &input, EntapDataPtrs &entap_data)
 *
 * Description           - Initializes member variables of Frame Selection
 *                       - Sets software type
 *
 * Notes                 - Called from EntapExecute, entry to Frame Selection
 *
 * @param input          - Input transcriptome (may have been from expression analysis)
 * @param entap_data     - Pointers to data needed during frame selection
 *
 * @return               - FrameSelection object
 *
 * =====================================================================
 */
FrameSelection::FrameSelection(std::string &input, EntapDataPtrs &entap_data) {
    FS_dprint("Spawn Object - FrameSelection");

    mQueryData        = entap_data.mpQueryData;
    mInPath           = input;
    mpFileSystem      = entap_data.mpFileSystem;
    mpUserInput       = entap_data.mpUserInput;

    mEntapDataPtrs = entap_data;

    mOutpath         = mpFileSystem->get_root_path();
    mOverwrite       = mpUserInput->has_input(INPUT_FLAG_OVERWRITE);
    mSoftwareFlag    = mpUserInput->get_user_input<ent_input_uint_t >(INPUT_FLAG_FRAME_SELECTION);

    mModOutDir   = PATHS(mOutpath, FRAME_SELECTION_OUT_DIR);
}


/**
 * =======================================================================
 * Function std::string FrameSelection::execute(std::string input)
 *
 * Description           - Handles overall Frame Selection process
 *                       - Runs software, parses and adds information to existing
 *                         data
 *
 * Notes                 - None
 *
 * @param input          - Input transcriptome (may have been from expression analysis)
 *
 * @return               - Path to new, frame selected, transcriptome
 *
 * ======================================================================
 */
std::string FrameSelection::execute(std::string input) {

    std::string output;                         // Absolute path to frame selected transcrtiptome
    EntapModule::ModVerifyData verify_data;     // Contains execution callback information
    std::unique_ptr<AbstractFrame> ptr;         // Pointer to frame selection object used during process


    mInPath = input;
    if (mOverwrite) mpFileSystem->delete_dir(mModOutDir);
    mpFileSystem->create_dir(mModOutDir);
    try {
        ptr = spawn_object();
        verify_data = ptr->verify_files();
        if (!verify_data.files_exist) {
            ptr->execute();
        }
        output = ptr->get_final_faa();
        ptr->parse();
        // If successful, set flags
        ptr->set_success_flags();

        ptr.reset();

        return output;
    } catch (const ExceptionHandler &e) {
        ptr.reset();
        throw e;
    }
}


/**
 * ======================================================================
 * Function std::unique_ptr<AbstractFrame> FrameSelection::spawn_object()
 *
 * Description           - Spawns object for specified frame selection software
 *
 * Notes                 - Currently supports GeneMarkS-T and TransDecoder
 *
 *
 * @return               - Frame Selection module object
 *
 * =====================================================================
 */
std::unique_ptr<AbstractFrame> FrameSelection::spawn_object() {
    switch (mSoftwareFlag) {
        case FRAME_GENEMARK_ST: {
            return std::unique_ptr<AbstractFrame>(new ModGeneMarkST(
                    mModOutDir,
                    mInPath,
                    mEntapDataPtrs
            ));
        }
        case FRAME_TRANSDECODER: {
            return std::unique_ptr<AbstractFrame>(new ModTransdecoder(
                    mModOutDir,
                    mInPath,
                    mEntapDataPtrs
            ));
        }
        default: {
            return std::unique_ptr<AbstractFrame>(new ModGeneMarkST(
                    mModOutDir,
                    mInPath,
                    mEntapDataPtrs
            ));
        }
    }
}

