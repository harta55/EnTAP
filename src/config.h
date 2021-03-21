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

#ifndef ENTAP_CONFIG_H
#define ENTAP_CONFIG_H

// Compile with boost libraries?
#ifndef USE_BOOST
//#define USE_BOOST   1
#endif

// Use EggNOG mapper (not supported, leaving for now)
#ifndef EGGNOG_MAPPER
//#define EGGNOG_MAPPER 1
#endif

// Compile with CURL? Will use wget command otherwise (not supported yet)
#ifndef USE_CURL
//#define USE_CURL    1
#endif

// Compile with using the Fast CSV Parser (required now)
#ifndef USE_FAST_CSV
#define USE_FAST_CSV  1
#endif

// Compile with ZLIB? Will use tar command otherwise (not supported yet)
#ifndef USE_ZLIB
//#define USE_ZLIB    1
#endif

//  Compile this in to run Unit tests
//#define UNIT_TESTS 1

// Comment this out if it is debug code
#define RELEASE_BUILD


#endif //ENTAP_CONFIG_H
