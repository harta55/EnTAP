/*
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

#ifndef ENTAP_VERSION_H
#define ENTAP_VERSION_H

#include "config.h"

#define LICENSE_YEAR_START      "2017"
#define LICENSE_YEAR_END        "2021"

// When changing Version ensure the EnTAP Database version/FTP is up-to-date (EntapDatabase.h)
#define MAJOR_VERSION     0
#define MINOR_VERSION     10
#define BUILD_VERSION     8

#define TO_STR2(x)             #x
#define TO_STR(x)              TO_STR2(x)

#ifndef RELEASE_BUILD
    #define ENTAP_VERSION_STR       (TO_STR(MAJOR_VERSION) "." TO_STR(MINOR_VERSION) "." \
                                    TO_STR(BUILD_VERSION) "-DEBUG")
#else
    #define ENTAP_VERSION_STR      (TO_STR(MAJOR_VERSION) "." TO_STR(MINOR_VERSION) "." \
                                    TO_STR(BUILD_VERSION))
#endif

#endif //ENTAP_VERSION_H
