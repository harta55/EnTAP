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

#ifndef ENTAP_ONTOLOGY_H
#define ENTAP_ONTOLOGY_H

#include "common.h"
#include "QueryData.h"
#include "EntapGlobals.h"
#include "QuerySequence.h"
#include "database/SQLDatabaseHelper.h"
#include "GraphingManager.h"
#include "ontology/AbstractOntology.h"
#include "EntapModule.h"

class AbstractOntology;

class Ontology {

public:

    void execute();
    Ontology(std::string, EntapDataPtrs &);

private:

    const std::string ONTOLOGY_OUT_PATH     = "ontology/";
    const std::string FINAL_ANNOT_FILE      = "final_annotations";
    const std::string FINAL_ANNOT_FILE_CONTAM = "final_annotations_contam";
    const std::string FINAL_ANNOT_FILE_NO_CONTAM = "final_annotations_no_contam";

    std::vector<ENTAP_HEADERS>      mEntapHeaders;
    std::vector<uint16>             mGoLevels;
    std::vector<uint16>             mSoftwareFlags;
    bool                            mIsOverwrite;
    std::string                     mOutpath;
    std::string                     mNewInput;
    std::string                     mOntologyDir;
    std::string                     mFinalOutputDir;
    QueryData                       *mpQueryData;
    FileSystem                      *mpFileSystem;
    UserInput                       *mpUserInput;
    EntapDataPtrs                   mEntapDataPtrs;
    std::vector<FileSystem::ENT_FILE_TYPES> mAlignmentFileTypes;

    void print_eggnog();
    std::unique_ptr<AbstractOntology> spawn_object(uint16&);
};


#endif //ENTAP_ONTOLOGY_H
