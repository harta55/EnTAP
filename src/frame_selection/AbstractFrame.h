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


#ifndef ENTAP_ABSTRACTFRAME_H
#define ENTAP_ABSTRACTFRAME_H

//*********************** Includes *****************************
#include "../QueryData.h"
#include "../EntapModule.h"
#include "../GraphingManager.h"
//**************************************************************

/**
 * ======================================================================
 * Class AbstractFrame
 *
 * Description          - Abstract class/functions for Frame Selection
 *                        portion of pipeline
 *                      - New modules will inherit from this class and
 *                        implement it's functions to fit into pipeline
 *
 * Notes                - None
 *
 *
 * @return              - None
 * ======================================================================
 */
class AbstractFrame : public EntapModule {

public:

    AbstractFrame(std::string &execution_stage_path, std::string &in_hits,
                  EntapDataPtrs &entap_data, std::string module_name,
                  std::vector<ENTAP_HEADERS> &module_headers);

    virtual ~AbstractFrame() = default;
    virtual ModVerifyData verify_files()=0;
    virtual void execute() = 0;
    virtual void parse() = 0;
    virtual bool set_version() = 0;

    virtual void set_success_flags() override ;

    std::string get_final_faa();


protected:
    const std::string FRAME_SELECTION_FIVE_FLAG     = "Partial 5 Prime";
    const std::string FRAME_SELECTION_THREE_FLAG    = "Partial 3 Prime";
    const std::string FRAME_SELECTION_COMPLETE_FLAG = "Complete";
    const std::string FRAME_SELECTION_INTERNAL_FLAG = "Internal";
    const std::string FRAME_SELECTION_LOST_FLAG     = "Lost";
    const std::string FRAME_SELECTION_FILENAME_PARTIAL   = "partial_genes";
    const std::string FRAME_SELECTION_FILENAME_COMPLETE  = "complete_genes";
    const std::string FRAME_SELECTION_FILENAME_INTERNAL  = "internal_genes";
    const std::string FRAME_SELECTION_FILENAME_LOST      = "sequences_removed";

    const std::string GRAPH_TITLE_FRAME_RESULTS     = "Frame_Selection_ORFs";
    const std::string GRAPH_FILE_FRAME_RESUTS       = "frame_results_pie.png";
    const std::string GRAPH_TEXT_FRAME_RESUTS       = "frame_results_pie.txt";
    const std::string GRAPH_TITLE_REF_COMPAR        = "Frame_Selected_Sequences";
    const std::string GRAPH_FILE_REF_COMPAR         = "removed_comparison_box.png";
    const std::string GRAPH_TEXT_REF_COMPAR         = "removed_comparison_box.txt";
    const std::string GRAPH_REJECTED_FLAG           = "Removed";
    const std::string GRAPH_KEPT_FLAG               = "Selected";

    const uint16 MINIMUM_KEPT_SEQUENCES = 0;        // Minimum required kept sequences before continuing pipeline

    std::string mFinalFaaPath;      // Absolute path to FAA (protein) file produced from Frame Selection process


    void frame_calculate_statistics();
};


#endif //ENTAP_ABSTRACTFRAME_H
