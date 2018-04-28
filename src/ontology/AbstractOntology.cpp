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


#include <boost/serialization/map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "../EntapGlobals.h"
#include "AbstractOntology.h"
#include "../ExceptionHandler.h"
#include "../database/EntapDatabase.h"


go_format_t AbstractOntology::parse_go_list(std::string list, EntapDatabase* database,char delim) {

    go_format_t output;
    std::string temp;
    std::vector<std::vector<std::string>>results;

    if (list.empty()) return output;
    std::istringstream ss(list);
    while (std::getline(ss,temp,delim)) {
        GoEntry term_info =
                database->get_go_entry(temp);
        output[term_info.category].push_back(temp + "-" + term_info.term +
                                             "(L=" + term_info.level + ")");
    }
    return output;
}