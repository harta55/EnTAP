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

#ifndef ENTAP_BASEREADER_H
#define ENTAP_BASEREADER_H


#include "../FileSystem.h"

class BaseReader {

protected:
    typedef enum {

        ERR_READ_OK=0,
        ERR_READ_FILE_OPEN,  // cannot open file
        ERR_READ_FILE_EMPTY, // file empty
        ERR_READ_FILE_PARSE,

    } READER_ERR;

    std::string print_err();
    BaseReader(std::string& file_path, FileSystem* fileSystem, FileSystem::ENT_FILE_TYPES file_type);
    virtual ~BaseReader(){};

    FileSystem                *_pFileSystem;
    std::string                _file_path;
    READER_ERR                 _reader_err;
    FileSystem::ENT_FILE_TYPES _file_type;
};


#endif //ENTAP_BASEREADER_H
