/******************************************************************
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2023, Alexander Hart, Dr. Jill Wegrzyn
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
 *******************************************************************/

#include "GFF_File.h"

GFF_File::GFF_File(FileSystem *fileSystem, QueryData *queryData, std::string &file_path) {
    mpFileSystem = fileSystem;
    mpQueryData = queryData;
    mGFFPath = file_path;
}

bool GFF_File::process_gff() {
    std::string line;
    std::string transcript_id;
    QuerySequence *current_query_sequence= nullptr;
    QuerySequence *previous_query_sequence= nullptr;
    bool ret = true;

    if (!mpFileSystem->file_exists(mGFFPath)) {
        mErrMessage = "ERROR GFF file does not exist at path: " + mGFFPath;
        return false;   // EXIT function
    }

    /*
     * Typical GFF file is tab-delim, example below:
     * scaffold_0001	AUGUSTUS	mRNA	15714	19241	0.61	-	.	ID=Pa1_00001.1;Parent=Pa1_00001
     *  we are interested in the 'mRNA' or 'transcript' lines of the GFF file, this will give us the ID
     *  after we have the correct line, we want to pull the start/stop information and the ID 9from 'ID=Pa1_0001.1'
     *
     *  GFF file can vary in format quite a bit, may not always have the same number of columns, so will try to parse it
     *      somewhat generically
     *
     * */

    std::ifstream in_file(mGFFPath);
    while (getline(in_file, line)) {
        if (line.empty()) continue;

        // Check if the line contains 'mRNA' or 'transcript', and skip if neither
        if ((line.find(TRANSCRIPT_ID_TAG_1) == std::string::npos) && (line.find(TRANSCRIPT_ID_TAG_2) == std::string::npos)) {
            continue;   // CONTINUE if this is not a line that contains the 'mRNA' or 'transcript' tag
        }

        // Check if we have known tags to determine transcript ID
        auto ind1 = line.find(TRANSCRIPT_ID_START);
        auto ind2 = line.find(TRANSCRIPT_ID_END);

        if (ind1 == std::string::npos || ind2 == std::string::npos) {
            continue;   // CONTINUE if we cannot find what we need to determine transcript ID
        }

        transcript_id = line.substr(ind1 + TRANSCRIPT_ID_START.length(),
                                    ind2-(ind1+TRANSCRIPT_ID_START.length()));

        current_query_sequence = mpQueryData->get_sequence(transcript_id);
        if (current_query_sequence != nullptr) {
            if (previous_query_sequence != nullptr) {
                previous_query_sequence->setMpDownstreamSequence(current_query_sequence);
                current_query_sequence->setMpUpstreamSequence(previous_query_sequence);
            }
            previous_query_sequence = current_query_sequence;
        } else {
            mErrMessage = "Unable to find sequence in the GFF file: " + transcript_id;
            ret = false;
            break;
        }
    }
    return ret;
}

const std::string &GFF_File::getMErrMessage() const {
    return mErrMessage;
}
