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

#include "NCBIEntrez.h"

const std::string NCBIEntrez::NCBI_DATABASE_TAXONOMY = "taxonomy";
const std::string NCBIEntrez::NCBI_ENTREZ_RETMODE_TEXT = "text";
const std::string NCBIEntrez::NCBI_ENTREZ_RETTYPE_XML  = "xml";

NCBIEntrez::NCBIEntrez(FileSystem *fileSystem) {
    mpFileSystem = fileSystem;

}

NCBIEntrez::~NCBIEntrez() {

}

bool NCBIEntrez::entrez_has_hits(NCBIEntrez::EntrezInput &entrezInput) {
    bool ret = false;
    std::string query;
    std::string output_file;
    std::string line;
    uint64 ind1;
    uint64 ind2;
    EntrezResults entrezResults;

    entrezInput.rettype = NCBI_ENTREZ_RETTYPE_COUNT;
    query = generate_query(entrezInput, ENTREZ_SEARCH);

    ret = mpFileSystem->download_url_data(query, output_file);

    // if successful download, parse
    if (ret) {
        if (mpFileSystem->file_exists(output_file)) {
            std::ifstream infile(output_file);
            while (std::getline(infile, line)) {
                // Look for count
                if (line.find(NCBI_ENTREZ_COUNT_XML_START) != std::string::npos) {
                    ind1 = line.find(NCBI_ENTREZ_COUNT_XML_START);
                    ind2 = line.find(NCBI_ENTREZ_COUNT_XML_END);
                    entrezResults.count = line.substr(ind1 + NCBI_ENTREZ_COUNT_XML_START.size(),ind2 - ind1 - NCBI_ENTREZ_COUNT_XML_START.size());
                    break;  // BREAK out of loop once we find count
                }
            }
            infile.close();
        } else {
            ret = false;
        }
    }
    return entrezResults.count != NCBI_ENTREZ_COUNT_ZERO;
}

bool NCBIEntrez::entrez_search(NCBIEntrez::EntrezInput &entrezInput, NCBIEntrez::EntrezResults &entrezResults) {
    bool ret=false;
    std::string output_file;
    std::string line;
    std::string temp;
    std::string query = generate_query(entrezInput, ENTREZ_SEARCH);
    uint64 ind1;
    uint64 ind2;

    // Download data
    ret = mpFileSystem->download_url_data(query, output_file);

    // if successful download, parse
    if (ret) {
        if (mpFileSystem->file_exists(output_file)) {
            std::ifstream infile(output_file);
            while (std::getline(infile, line)) {

                // Look for ID
                if (line.find(NCBI_ENTREZ_UID_XML_START) != std::string::npos) {
                    ind1 = line.find(NCBI_ENTREZ_UID_XML_START);
                    ind2 = line.find(NCBI_ENTREZ_UID_XML_END);
                    temp = line.substr(ind1 + NCBI_ENTREZ_UID_XML_START.size(),ind2 - ind1 - NCBI_ENTREZ_UID_XML_START.size());
                    entrezResults.uid_list.push_back(temp);
                }

                // Look for count
                if (line.find(NCBI_ENTREZ_COUNT_XML_START) != std::string::npos) {
                    ind1 = line.find(NCBI_ENTREZ_COUNT_XML_START);
                    ind2 = line.find(NCBI_ENTREZ_COUNT_XML_END);
                    entrezResults.count = line.substr(ind1 + NCBI_ENTREZ_COUNT_XML_START.size(),ind2 - ind1 - NCBI_ENTREZ_COUNT_XML_START.size());
                }
            }
            infile.close();
        } else {
            ret = false;
        }
    }
    return ret;
}

std::string NCBIEntrez::generate_query(NCBIEntrez::EntrezInput &entrezInput, ENTREZ_TYPES type) {

    std::string final_url;

    switch (type) {

        case ENTREZ_SEARCH:
            final_url = ESEARCH_BASE_URL;
            break;

        case ENTREZ_SUMMARY:
            final_url = ESUMMARY_BASE_URL;
            break;

        case ENTREZ_FETCH:
            final_url = EFETCH_BASE_URL;
            break;
    }

    if (!entrezInput.term.empty()) {
        process_term(entrezInput.term); // NO spaces in term
        final_url += "&" + NCBI_ENTREZ_TERM + entrezInput.term;
    }

    if (!entrezInput.database.empty()) {
        final_url += "&" + NCBI_ENTREZ_DATABASE + entrezInput.database;
    }

    if (!entrezInput.retmode.empty()) {
        final_url += "&" + NCBI_ENTREZ_RETMODE + entrezInput.retmode;
    }

    if (!entrezInput.rettype.empty()) {
        final_url += "&" + NCBI_ENTREZ_RETTYPE + entrezInput.rettype;
    }

    if (!entrezInput.uid_list.empty()) {
        final_url += "&" + NCBI_ENTREZ_UID;

        for (auto &uid : entrezInput.uid_list) {
            if (!uid.empty()) {
                final_url += uid + ",";
            }
        }
        final_url.pop_back(); // remove trailing ','
    }
    return final_url;
}

void NCBIEntrez::process_term(std::string &term) {
    STR_REPLACE(term, ' ', '_');
}
