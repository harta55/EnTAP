/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
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
#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <boost/regex.hpp>
#include <iomanip>
#include "FrameSelection.h"
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapConfig.h"
#include "EntapGlobals.h"
#include "frame_selection/ModGeneMarkST.h"
#include "QueryData.h"
//**************************************************************


/**
 * ======================================================================
 * Function FrameSelection(std::string &input, std::string &out,
 *                             boost::program_options::variables_map &user_flags,
 *                             GraphingManager *graphingManager)
 *
 * Description           - Initializes member variables of Frame Selection
 *                       - Sets software type
 *
 * Notes                 - Called from EntapExecute
 *
 * @param input          - Input transcriptome (may have been from expression analysis)
 * @param out            - EnTAP main output directory
 * @param user_flag      - Boost map of user inputs
 * @param graphingManager- Pointer to graphing manager
 *
 * @return               - FrameSelection object
 *
 * =====================================================================
 */
FrameSelection::FrameSelection(std::string &input, std::string &out,
                               boost::program_options::variables_map &user_flags,
                               GraphingManager *graphingManager, QueryData *QUERY_DATA) {
    print_debug("Spawn object - FrameSelection");
    _graphingManager = graphingManager;
    _QUERY_DATA = QUERY_DATA;
    _exe_path        = GENEMARK_EXE;
    _outpath         = out;
    _inpath          = input;
    _overwrite       = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _software_flag   = 0;
    _processed_path  = (boostFS::path(out) / boostFS::path(FRAME_SELECTION_OUT_DIR) /
                        boostFS::path(PROCESSED_DIR)).string();
    _figure_path     = (boostFS::path(out) / boostFS::path(FRAME_SELECTION_OUT_DIR) /
                        boostFS::path(FIGURE_DIR)).string();
    _frame_outpath   = (boostFS::path(out) / boostFS::path(FRAME_SELECTION_OUT_DIR)).string();
    SOFTWARE = static_cast<FrameSoftware >(_software_flag);
}


/**
 * ======================================================================
 * Function std::string FrameSelection::execute(std::string input,
 *                                              std::map<std::string,QuerySequence> &SEQUENCES)
 *
 * Description           - Handles overall Frame Selection process
 *                       - Runs software, parses and adds information to existing
 *                         data
 *
 * Notes                 - None
 *
 * @param input          - Input transcriptome (may have been from expression analysis)
 * @param SEQUENCES      - Sequence map of information thus far
 *
 * @return               - Path to new transcriptome
 *
 * =====================================================================
 */
std::string FrameSelection::execute(std::string input) {

    std::string output;
    std::pair<bool, std::string> verify_pair;

    _inpath = input;
    if (_overwrite) boostFS::remove_all(_frame_outpath);
    boostFS::create_directories(_frame_outpath);
    boostFS::create_directories(_processed_path);
    boostFS::create_directories(_figure_path);
    try {
        std::unique_ptr<AbstractFrame> ptr = spawn_object();
        verify_pair = ptr->verify_files();
        if (!verify_pair.first) {
            output = ptr->execute();
        } else output = verify_pair.second;
        ptr->parse();
        ptr.release();
        return output;
    } catch (const ExceptionHandler &e) {throw e;}
}


/**
 * ======================================================================
 * Function std::unique_ptr<AbstractFrame> FrameSelection::spawn_object()
 *
 * Description           - Spawns object for specified frame selection software
 *
 * Notes                 - Currently only GeneMarkS-T support
 *
 *
 * @return               - Frame Selection module object
 *
 * =====================================================================
 */
std::unique_ptr<AbstractFrame> FrameSelection::spawn_object() {
    // Handle any special conditions for each software

    switch (SOFTWARE) {
        case GENEMARKST:
            return std::unique_ptr<AbstractFrame>(new ModGeneMarkST(
                    _exe_path, _outpath, _inpath, _processed_path, _figure_path,
                    _frame_outpath, _graphingManager, _QUERY_DATA
            ));
        default:
            return std::unique_ptr<AbstractFrame>(new ModGeneMarkST(
                    _exe_path, _outpath, _inpath, _processed_path, _figure_path,
                    _frame_outpath, _graphingManager, _QUERY_DATA
            ));
    }
}

