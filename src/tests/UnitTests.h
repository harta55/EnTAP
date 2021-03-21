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

#ifndef ENTAP_UNITTESTS_H
#define ENTAP_UNITTESTS_H

#include "../config.h"

#ifdef UNIT_TESTS
#include "../FileSystem.h"
#include "../database/EntapDatabase.h"

class UnitTests {

public:
    UnitTests();
    ~UnitTests();

    void execute_tests();

protected:
    FileSystem* mpFileSystem;
    std::string mRootDirectory;

    void TestEntapDatabase();
    void UTEntapDatabase_00(EntapDatabase *entapDatabase, EntapDatabase::DATABASE_TYPE type);
    void UTEntapDatabase_01(EntapDatabase *entapDatabase, EntapDatabase::DATABASE_TYPE type);
    void UTEntapDatabase_02(EntapDatabase *entapDatabase, EntapDatabase::DATABASE_TYPE type);
    void UTEntapDatabase_03(EntapDatabase *entapDatabase, EntapDatabase::DATABASE_TYPE type);

    void TestQueryData();

    void TestEggnogDatabase();

};

#endif //ENTAP_UNITTESTS_H
#endif