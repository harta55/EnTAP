/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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

//*********************** Common Includes ******************************
#include <vector>
#include <string>
#include <list>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <cmath>
#include <exception>
#include <queue>
#include <thread>
#include <unordered_map>
#include <chrono>
#include <ios>
#include <iterator>
//**************************************************************


//******************* Defines/Macros ***************************
#define LOWERCASE(x)        std::transform(x.begin(), x.end(), x.begin(), ::tolower)
#define STR_COUNT(x,y)      std::count(x.begin(), x.end(), y)
#define STR_REPLACE(x,y,z)  std::replace(x.begin(), x.end(), y, z)
#define STR_ERASE(x,y)      x.erase(std::remove(x.begin(), x.end(), y),x.end())
#define FIND_VECT(x,y)      std::find(y.begin(), y.end(), x) != y.end() // find x in y
#define SAFE_DELETE(x)      if (x) delete x; x = nullptr
#define ENTAP_PERCENT       (100.0f)
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

typedef std::pair<std::string,std::string>    pair_str_t;
typedef std::vector<std::string>              vect_str_t;
typedef std::vector<std::vector<std::string>> vect_vect_str_t;
typedef std::vector<uint16>                   vect_uint16_t;
typedef std::vector<fp64>                     vect_fp64_t;
typedef std::set<std::string>                 set_str_t;
//**************************************************************

#endif //ENTAP_COMMON_H
