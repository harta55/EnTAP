/*
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
*/

#ifndef ENTAP_MODHORIZONTALGENETRANSFERDIAMOND_H
#define ENTAP_MODHORIZONTALGENETRANSFERDIAMOND_H


#include "AbstractHorizontalGeneTransfer.h"
#include "../common.h"
#include "../similarity_search/AbstractSimilaritySearch.h"
#include "../similarity_search/ModDiamond.h"

class ModHorizontalGeneTransferDiamond : public AbstractHorizontalGeneTransfer{

public:
    ModHorizontalGeneTransferDiamond(std::string &execution_stage_path, std::string &fasta_path, EntapDataPtrs &entap_data);
    ~ModHorizontalGeneTransferDiamond()=default;

    // ModEntap overrides
    virtual ModVerifyData verify_files() override;
    virtual void execute() override ;
    virtual void parse() override ;
    static bool is_executable(std::string& exe);
    virtual bool set_version() override;

private:
    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;

    std::unique_ptr<AbstractSimilaritySearch> mModDiamondDonor;
    std::unique_ptr<AbstractSimilaritySearch> mModDiamondRecipient;
};


#endif //ENTAP_MODHORIZONTALGENETRANSFERDIAMOND_H
