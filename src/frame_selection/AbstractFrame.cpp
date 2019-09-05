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

#include "AbstractFrame.h"

/**
 * ======================================================================
 * Function AbstractFrame(std::string &execution_stage_path, std::string &in_hits,
 *                        EntapDataPtrs &entap_data, std::string module_name,
                          std::string &exe)
 *
 * Description          - Constructor for Abstract frame selection class
 *                      - Initializes protected member variables for
 *                        expression modules
 *
 * Notes                - Constructor
 *
 * @param execution_stage_path - Absolute path to the output directory for this stage (Frame Selection)
 * @param in_hits              - Absolute path to input transcriptome
 * @param entap_data           - Pointers to necessary entap data for frame selection
 * @param module_name          - Name of this software module
 * @param exe                  - Execution method (i.e. executable)
 *
 * @return              - AbstractFrame object
 * ======================================================================
 */
AbstractFrame::AbstractFrame(std::string &execution_stage_path, std::string &in_hits,
                             EntapDataPtrs &entap_data, std::string module_name, std::string &exe)
: EntapModule(execution_stage_path, in_hits, entap_data, module_name, exe) {

    mExecutionState = FRAME_SELECTION;
}