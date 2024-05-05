/******************************************************************
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2024, Alexander Hart, Dr. Jill Wegrzyn
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
 *******************************************************************/

#ifndef ENTAP_GFF_FILE_H
#define ENTAP_GFF_FILE_H


#include "../FileSystem.h"
#include "../QueryData.h"

class GFF_File {
public:
    GFF_File(FileSystem *fileSystem, QueryData *queryData, std::string &file_path);
    bool process_gff();

private:
    FileSystem *mpFileSystem;   // Pointer to EnTAP filesystem
    QueryData *mpQueryData;
    std::string mGFFPath;       // File path to GFF File
    std::string mErrMessage;
public:
    const std::string &getMErrMessage() const;

private:

    const uint8 GFF_COL_NUMBER = 9;    // 'Standard' GFF file column number
    const std::string TRANSCRIPT_ID_TAG_1 = "mRNA";
    const std::string TRANSCRIPT_ID_TAG_2 = "transcript";
    const std::string TRANSCRIPT_ID_START = "ID=";
    const std::string TRANSCRIPT_ID_END   = ";";

};


#endif //ENTAP_GFF_FILE_H
