/*
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

#include "BaseReader.h"

std::string BaseReader::print_err() {
    std::string err_out;

    switch (_reader_err) {
        case ERR_READ_OK:
            err_out = "No error";
            break;

        case ERR_READ_FILE_OPEN:
            err_out = "Unable to open file";
            break;

        case ERR_READ_FILE_EMPTY:
            err_out = "File empty";
            break;

        default:
            err_out = "Unknown error code";
            break;
    }
    return err_out;
}

BaseReader::BaseReader(std::string &file_path,
                       FileSystem *fileSystem,
                       FileSystem::ENT_FILE_TYPES file_type) {
    _file_path = file_path;
    _pFileSystem = fileSystem;
    _file_type = file_type;
}
