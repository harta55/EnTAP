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

#include "NCBIEntrez.h"
#include <regex>
#include "NCBIDatabase.h"

const std::string NCBIEntrez::NCBI_DATABASE_TAXONOMY = "taxonomy";
const std::string NCBIEntrez::NCBI_DATABASE_PROTEIN = "protein";

const std::string NCBIEntrez::NCBI_ENTREZ_RETMODE_TEXT = "text";

const std::string NCBIEntrez::NCBI_ENTREZ_RETTYPE_XML  = "xml";
const std::string NCBIEntrez::NCBI_ENTREZ_RETTYPE_GP  = "gp";

NCBIEntrez::NCBIEntrez(FileSystem *fileSystem) {
    mpFileSystem = fileSystem;
}

NCBIEntrez::~NCBIEntrez() {

}

bool NCBIEntrez::entrez_has_hits(EntrezInput &entrezInput) {
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

bool NCBIEntrez::entrez_fetch(EntrezInput& entrezInput, EntrezResults& entrezResults) {
    std::string output_file;
    std::string line;           // line read from file downloaded from NCBI

    if (entrezInput.uid_list.empty()) return false;
    std::string query = generate_query(entrezInput, ENTREZ_FETCH);
    if (query.empty()) return false;

    // TODO will need to include CURL in project, downloading to file now to test
    // Download data
    bool ret = mpFileSystem->download_url_data(query, output_file);
    if (ret && mpFileSystem->file_exists(output_file))
    {
        // File successfully downloaded, parse
        if (entrezInput.rettype == NCBI_ENTREZ_RETTYPE_GP) {
            return parse_ncbi_gp_file(entrezInput, entrezResults, output_file);
        } else {
            return false;
        }
    } else {
        return false;
    }
}

std::string NCBIEntrez::generate_query(EntrezInput &entrezInput, ENTREZ_TYPES type) {

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

/*
 * GP Flat Files from NCBI generally start with the following format:
 *
 * (+)
 * DEFINITION  carbonic anhydrase 2 [Cimex lectularius].
 * ACCESSION   XP_014245616
 * VERSION     XP_014245616.1
 * DBLINK      BioProject: PRJNA298750
 * DBSOURCE    REFSEQ: accession XM_014390130.2
 * KEYWORDS    RefSeq.
 * SOURCE      Cimex lectularius (bed bug)
 * ........................ (a bit further down in file)
 *      CDS             1..298
 *                    /gene="LOC106664428"
 *                    /coded_by="XM_014390130.2:154..1050"
 *                    /db_xref="GeneID:106664428"
 *
 *  Unfortunately access to the database requires manual parsing to get specific data
 */
bool NCBIEntrez::parse_ncbi_gp_file(EntrezInput& entrezInput, EntrezResults& entrezResults, std::string &output_file) {
    std::string line;
    std::smatch match;
    std::string current_sequence;
    EntrezEntryData entrez_entry_data;
    std::ifstream infile(output_file);
    const std::string TEST_REGEX = R"(LOCUS\s*(\S+)\s)";
    const std::string TEXT_REGEX_GENE = "\\s+\\/db_xref=\"GeneID:(.+)\"";

    // NOTE with NCBI data versions, depending on what database version we are searching against,
    //  it could be inconsistent with the latest version on NCBI. These versions are denoted after the '.'
    //  example: XP_014245616.1 is verison 1. Theoretically it could be XP_014245616.2 at some point. We
    //  also want to preserve the version that the user searched against.
    // Because of this, we are going to map  the Users data to what we find in NCBI, extremely inefficient...
    std::unordered_map<std::string, std::string> ncbi_id_mappings;
    for (std::string& entry : entrezInput.uid_list) {
        if (entry.empty()) continue;
        std::string reformatted = entry.substr(0,entry.find('.'));
        ncbi_id_mappings.emplace(reformatted, entry);
    }

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        // Find what sequence we are on (it is possible to pull data for multiple at once
        if (current_sequence.empty()) {
            if (std::regex_search(line, match, std::regex(TEST_REGEX))) {
                current_sequence = std::string(match[1]);
                entrez_entry_data = {};
            }
        } else {
            // Check if this line contains any data we want
            for (ENTREZ_DATA_TYPES data_type : entrezInput.data_types) {
                switch (data_type) {
                    case ENTREZ_DATA_GENEID:
                        if (std::regex_search(line, match, std::regex(TEXT_REGEX_GENE))) {
                            entrez_entry_data.geneid = std::string(match[1]);
                        }
                        break;
                    default:
                        break;

                }
            }

            // Check if have finished with this sequences
            //  NCBI denotes this as '//'
            if (line == "//") {
                entrezResults.entrez_results.emplace(ncbi_id_mappings.at(current_sequence), entrez_entry_data);
                current_sequence = "";
            }
        }
    }
    infile.close();
    return true;
}
