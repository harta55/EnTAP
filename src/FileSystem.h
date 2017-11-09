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

#ifndef ENTAP_FILESYSTEM_H
#define ENTAP_FILESYSTEM_H


//*********************** Includes *****************************
#include "common.h"
//**************************************************************


//***************** Global Prototype Functions *****************
#if 0
void FS_open_out(std::string &, std::ofstream &);
#endif
bool FS_file_is_open(std::ofstream&);
void FS_close_file(std::ofstream&);
void FS_dprint(std::string);
void FS_print_stats(std::string &msg);
bool FS_file_test_open(std::string&);
bool FS_file_exists(std::string);
bool FS_file_empty(std::string);
bool FS__delete_file(std::string);
bool FS_directory_iterate(bool, std::string&);
//**************************************************************

const std::string TXT_EXT = ".txt";


#endif //ENTAP_FILESYSTEM_H
