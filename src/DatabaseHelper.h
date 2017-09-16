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


#ifndef ENTAP_DATABASEHELPER_H
#define ENTAP_DATABASEHELPER_H

#include <iostream>
#include <vector>
#include "sqlite3.h"


class DatabaseHelper {

typedef std::vector<std::vector<std::string>> query_struct;

public:
    DatabaseHelper();
    ~DatabaseHelper();
    bool open(std::string file);
    void close();
    query_struct query(char* query);

private:
    sqlite3 *_database;
};



#endif //ENTAP_DATABASEHELPER_H
