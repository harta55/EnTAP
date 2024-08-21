/*
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
*/

#ifndef NCBIDATABASE_H
#define NCBIDATABASE_H
#include "../QueryAlignment.h"


class NCBIDatabase {

public:
 ~NCBIDatabase();
 NCBIDatabase(FileSystem *pfile_system);

 void get_ncbi_data(SimSearchAlignment* alignment);
private:
 bool mContinueProcessing;
 bool mProcessRemaining;
 std::list<std::thread::id> mRunningThreads;
 std::queue<SimSearchAlignment*> mProcessingQueue;
 NCBIEntrez *mpNCBIEntrez;

 void main_loop();
 void get_entrez_ncbi_data(const std::string& query_ids, std::unordered_map<std::string, SimSearchAlignment*> bulk_query_map);
 uint16 ENTREZ_QUEUE_SIZE = 100; // Bulk processing size for entrez data
};



#endif //NCBIDATABASE_H
