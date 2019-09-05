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


#ifndef ENTAP_ABSTRACTFRAME_H
#define ENTAP_ABSTRACTFRAME_H

//*********************** Includes *****************************
#include "../QueryData.h"
#include "../EntapModule.h"
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
                  EntapDataPtrs &entap_data, std::string module_name, std::string &exe);

    virtual ~AbstractFrame() = default;
    virtual ModVerifyData verify_files()=0;
    virtual void execute() = 0;
    virtual void parse() = 0;
    virtual std::string get_final_faa() = 0;

protected:
    const uint16 MINIMUM_KEPT_SEQUENCES = 0;        // Minimum required kept sequences before continuing pipeline

};


#endif //ENTAP_ABSTRACTFRAME_H
