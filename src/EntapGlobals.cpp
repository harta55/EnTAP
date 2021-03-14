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


//*********************** Includes *****************************
#include "EntapGlobals.h"
//**************************************************************

// For C++11 use std::to_string(arg)
std::string float_to_string(fp64 val) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << val;
    return ss.str();
}

std::string float_to_sci(fp64 val, int precision) {
    std::stringstream ss;
    ss << std::fixed <<
          std::setprecision(precision) <<
          std::scientific              <<
          val;
    return ss.str();
}

vect_str_t split_string(std::string sequences, char delim) {
    vect_str_t split_vals;

    // Remove newline if exists
    sequences.erase(std::remove(sequences.begin(), sequences.end(), '\n'), sequences.end());

    std::istringstream iss(sequences);
    std::string val;
    while(std::getline(iss, val, delim)) {
        split_vals.push_back(val);
    }
    return split_vals;
}

std::string get_cur_time() {
    std::chrono::time_point<std::chrono::system_clock> current;
    std::time_t time;

    current = std::chrono::system_clock::now();
    time = std::chrono::system_clock::to_time_t(current);
    std::string out_time(std::ctime(&time));
    return out_time.substr(0,out_time.length()-1);
}

std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}