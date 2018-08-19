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

#if 0

#ifndef ENTAP_DELIMREADER_H
#define ENTAP_DELIMREADER_H


#include "BaseReader.h"
#include "../config.h"

class DelimReader : public BaseReader {

public:
    DelimReader(std::string& file_path, FileSystem *filesystem, FileSystem::ENT_FILE_TYPES,
        char delim);
    ~DelimReader();
    std::pair<bool, vect_vect_str_t> parse_file();

private:
    char _delim;
    bool _success;
};


#endif //ENTAP_DELIMREADER_H

#endif