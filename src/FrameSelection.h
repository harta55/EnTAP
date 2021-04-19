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


#ifndef ENTAP_FRAMESELECTION_H
#define ENTAP_FRAMESELECTION_H

//*********************** Includes *****************************
#include "frame_selection/AbstractFrame.h"
#include "QueryData.h"
#include "common.h"
//**************************************************************

class AbstractFrame;

// Supported frame selection software
enum FRAME_SELECTION_SOFTWARE {
    FRAME_UNUSED=0,
    FRAME_GENEMARK_ST,
    FRAME_TRANSDECODER,
    FRAME_SOFTWARE_COUNT
};

class FrameSelection {

public:
    std::string execute(std::string);
    FrameSelection(std::string&, EntapDataPtrs&);

private:

    const std::string FRAME_SELECTION_OUT_DIR       = "frame_selection/";   // Frame selection directory name

    std::string      mModOutDir;    // Absolute path to Frame Selection out directory
    std::string      mInPath;       // Absolute path to input transcriptome
    std::string      mOutpath;      // Absolute path to root output directory
    bool             mOverwrite;    // TRUE if old files should be overwritten
    uint16           mSoftwareFlag; // Type of software to use during Frame Selection
    QueryData        *mQueryData;           // Pointer to query data
    FileSystem       *mpFileSystem;         // Pointer to EnTAP filesystem
    UserInput        *mpUserInput;          // Pointer to user input flags
    EntapDataPtrs    mEntapDataPtrs;        // Contains all data pointers needed during frame selection

    std::unique_ptr<AbstractFrame> spawn_object();
};


#endif //ENTAP_FRAMESELECTION_H
