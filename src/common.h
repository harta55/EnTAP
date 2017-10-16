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

#ifndef ENTAP_COMMON_H
#define ENTAP_COMMON_H

//*********************** Includes ******************************
#include <vector>
#include <string>
#include <list>
#include <fstream>
#include <iostream>
//**************************************************************


//******************* Defines/Macros ***************************
#define LOWERCASE(x)        std::transform(x.begin(), x.end(), x.begin(), ::tolower)

//**************************************************************


//*********************** Typedefs *****************************

typedef unsigned char       uint8;
typedef signed char         int8;
typedef unsigned short      uint16;
typedef signed short        int16;
typedef unsigned int        uint32;
typedef signed int          int32;
typedef unsigned long long  uint64;
typedef signed long long    int64;
typedef float               fp32;
typedef double              fp64;
//**************************************************************

#endif //ENTAP_COMMON_H