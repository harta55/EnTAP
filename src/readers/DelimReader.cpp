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
#include <csv.h>
#include "DelimReader.h"

DelimReader::DelimReader(std::string &file_path, FileSystem *filesystem,
                         FileSystem::ENT_FILE_TYPES file_type, char delim) :
BaseReader(file_path, filesystem, file_type) {
    FS_dprint("Spawn Object - DelimReader");

    _delim = delim;
}

DelimReader::~DelimReader() {
    FS_dprint("Kill Object - DelimReader");

}

std::pair<bool, vect_vect_str_t> DelimReader::parse_file() {
    uint32 column_num=0;
    vect_vect_str_t out_data;
    bool valid_parse = true;

    // Ensure file is valid before continuing
    if (_pFileSystem->get_file_status(_file_path) != 0)
        return std::make_pair(false, out_data);

    std::ifstream infile(_file_path);

    // Get number of columns in file
    if (infile.good()) {
        std::string first_line;
        getline(infile, first_line);
        for (uint32 i ; i < first_line.size(); i++) {
            if (first_line.at(i) == _delim) column_num++;
        }
        column_num++;
    }
    // Reset back to beginning of file
    infile.clear();
    infile.seekg(0, std::ios::beg);

#ifdef USE_FAST_CSV
    FS_dprint("Using Fast CSV to parse file...");



#endif


    return std::make_pair(false, out_data);
}
#endif