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

#ifndef ENTAP_ENTAPEXECUTE_H
#define ENTAP_ENTAPEXECUTE_H

//*********************** Includes *****************************
#include "Ontology.h"
#include "QuerySequence.h"
#include "QueryData.h"
#include <list>
#include <boost/program_options/variables_map.hpp>
#include <map>
#include <unordered_map>
#include <queue>

//**************************************************************

class QueryData;
namespace entapExecute {

    // ********************** Global Constants *********************
    const std::string ENTAP_OUTPUT         = "entap_out/";

    //**************************************************************


    // ****************** Global Prototype Functions ***************
    void execute_main(UserInput*, FileSystem*);
    //**************************************************************


}


#endif //ENTAP_ENTAPEXECUTE_H
