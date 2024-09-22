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

#include "NCBIDatabase.h"
#include "../QueryAlignment.h"
#include <thread>
#include <unistd.h>


NCBIDatabase::~NCBIDatabase()
{
    delete mpNCBIEntrez;
}

NCBIDatabase::NCBIDatabase(FileSystem *pfile_system)
{
    mContinueProcessing = true;
    mProcessRemaining = false;
    mpNCBIEntrez = new NCBIEntrez(pfile_system);
    mNCBIAccessionType = NCBI_ACCESSION_ENTREZ;
    // std::thread thNCBIData(&NCBIDatabase::main_loop, this);
}

void NCBIDatabase::get_ncbi_data(SimSearchAlignment* alignment)
{
    mProcessingQueue.push(alignment);
}

NCBIDataResults_t NCBIDatabase::get_ncbi_data(const vect_str_t &ncbi_accession) const {
    if (ncbi_accession.empty()) return {};

    switch (mNCBIAccessionType) {

        case NCBI_ACCESSION_ENTREZ: {
            if (mpNCBIEntrez == nullptr) return {}; // RETURN empty data
            NCBIEntrez::EntrezResults entrez_results;
            NCBIEntrez::EntrezInput entrez_input;
            NCBIDataResults_t ncbi_data_results;
            entrez_input.database = NCBIEntrez::NCBI_DATABASE_PROTEIN;
            entrez_input.rettype = NCBIEntrez::NCBI_ENTREZ_RETTYPE_GP;
            entrez_input.data_types = {NCBIEntrez::ENTREZ_DATA_GENEID};
            entrez_input.uid_list = ncbi_accession;

            if (mpNCBIEntrez->entrez_fetch(entrez_input, entrez_results)) {
                // NCBI database pull was successful
                return entrez_results.entrez_results;
            }
            break;
        }

        case NCBI_ACCESSION_API:
        default:
            break;
    }
    return {};
}

void NCBIDatabase::main_loop()
{
    while(mContinueProcessing)
    {
        if ((mProcessingQueue.size() >= ENTREZ_QUEUE_SIZE) || (mProcessRemaining))
        {
            std::string query_ids;
            std::unordered_map<std::string, SimSearchAlignment*> bulk_query_map;
            for (uint16 i=0; i<=ENTREZ_QUEUE_SIZE; i++)
            {
                if (!mProcessingQueue.empty())
                {
                    query_ids += mProcessingQueue.front()->get_results()->sseqid + ',';
                    bulk_query_map.emplace(mProcessingQueue.front()->get_results()->sseqid,
                        mProcessingQueue.front());
                } else
                {
                    break;
                }
            }
            // Start thread to pull NCBI data
            query_ids.pop_back();// remove trailing ','
            std::thread th(&NCBIDatabase::get_entrez_ncbi_data,this, query_ids, bulk_query_map);
            mRunningThreads.push_back(th.get_id());
        }
        sleep(5);
    }
}

void NCBIDatabase::get_entrez_ncbi_data(const std::string& query_ids, std::unordered_map<std::string, SimSearchAlignment*> bulk_query_map)
{
    if (mpNCBIEntrez == nullptr) return;    // RETURN with invalid entrez pointer
    SimSearchAlignment *alignment;
    NCBIEntrez::EntrezResults entrez_results;
    NCBIEntrez::EntrezInput entrez_input;
    entrez_input.database = NCBIEntrez::NCBI_DATABASE_PROTEIN;
    entrez_input.rettype = NCBIEntrez::NCBI_ENTREZ_RETTYPE_GP;
    entrez_input.uid_list = split_string(query_ids, ',');   // Redo code here to prevent unecessary manipulation

    if (mpNCBIEntrez->entrez_fetch(entrez_input, entrez_results))
    {
        for (auto &results : entrez_results.entrez_results) {
            alignment = bulk_query_map.at(results.first);
            if (alignment != nullptr) {
                //alignment->get_results().
            }
        }
    }
}
