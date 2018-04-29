/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
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
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"
#include "GraphingManager.h"
#include "frame_selection/AbstractFrame.h"
#include "QueryData.h"
#include "common.h"
//**************************************************************

class AbstractFrame;
class FrameSelection {


public:
    std::string execute(std::string);
    FrameSelection(std::string&, EntapDataPtrs&);


private:

    const std::string FRAME_SELECTION_OUT_DIR       = "frame_selection/";

    std::string      _frame_outpath;
    std::string      _exe_path;
    std::string      _inpath;
    std::string      _outpath;
    bool             _overwrite;
    short            _software_flag;
    GraphingManager  *_pGraphingManager;
    QueryData        *_QUERY_DATA;
    FileSystem       *_pFileSystem;
    UserInput        *_pUserInput;
    EntapDataPtrs    _entap_data_ptrs;

    std::unique_ptr<AbstractFrame> spawn_object();
};


#endif //ENTAP_FRAMESELECTION_H
