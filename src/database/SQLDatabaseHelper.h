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


#ifndef ENTAP_DATABASEHELPER_H
#define ENTAP_DATABASEHELPER_H

#include <iostream>
#include "../common.h"
#include "sqlite3.h"


class SQLDatabaseHelper {

public:
    typedef std::vector<std::vector<std::string>> query_struct;

    SQLDatabaseHelper();
    ~SQLDatabaseHelper();
    bool open(std::string file);
    bool create(std::string file);
    bool execute_cmd(char*);
    void close();
    query_struct query(char* query);

    // change to template
    std::string format_container(std::set<std::string> &in_cont);
    std::string format_string(std::string& str, char delim);

private:
    sqlite3 *mpDatabase;
};



#endif //ENTAP_DATABASEHELPER_H
