/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2020, Alexander Hart, Dr. Jill Wegrzyn
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

#ifndef ENTAP_MODBUSCO_H
#define ENTAP_MODBUSCO_H


#include "AbstractOntology.h"

class ModBUSCO : public AbstractOntology {
public:
    ModBUSCO(std::string &in_hits, std::string &ont_out, EntapDataPtrs &entap_data);
    ~ModBUSCO();
    virtual void execute() override;
    virtual void parse() override;
    virtual ModVerifyData verify_files() override;
    static bool is_executable(std::string &exe);

private:
    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;

    // WARNING current as of version 3.0.2, subject to change!!
    const std::string BUSCO_INPUT_IN       = "--in";        // Specify input transcriptome
    const std::string BUSCO_INPUT_OUTPUT   = "-o";          // Specify output flag/directory (sort of works?)
    const std::string BUSCO_INPUT_RUN_TYPE = "-m";          // Specify execution type (prot, tran)
    const std::string BUSCO_INPUT_DATABASE = "-l";          // Specify path to BUSCO database
    const std::string BUSCO_INPUT_CPU      = "--cpu";       // Specify number of threads to BUSCO

    // BUSCO prepends this tag before all directory creations (-o command)
    const std::string BUSCO_PREPEND_TAG    = "run_";

    const std::string BUSCO_RUN_TYPE_TRAN  = "tran";
    const std::string BUSCO_RUN_TYPE_PROT  = "prot";

    std::string mExePath;      // Path to BUSCO executable
    std::string mRunType;      // Execution type for BUSCO (prot or tran based on blastp)
    vect_str_t  mDatabasePaths; // Vector of BUSCO databases we want to use
    fp64        mEval;         // Eval to use during BUSCO search


};


#endif //ENTAP_MODBUSCO_H
