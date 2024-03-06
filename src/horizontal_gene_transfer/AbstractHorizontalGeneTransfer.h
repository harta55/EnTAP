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

#ifndef ENTAP_ABSTRACTHORIZONTALGENETRANSFER_H
#define ENTAP_ABSTRACTHORIZONTALGENETRANSFER_H


#include "../EntapModule.h"
#include "../QueryData.h"

class AbstractHorizontalGeneTransfer : public EntapModule {

public:

    AbstractHorizontalGeneTransfer(std::string &execution_stage_path, std::string &in_hits,
                                   EntapDataPtrs &entap_data, std::string mod_name,
                                   std::vector<ENTAP_HEADERS> &module_headers);

    ~AbstractHorizontalGeneTransfer() = default;
    virtual ModVerifyData verify_files()=0;
    virtual void execute() = 0;
    virtual void parse() = 0;
    virtual void set_success_flags() override ;
    virtual bool set_version() = 0;

protected:
    std::string                     mInputLineage;
    std::string                     mInputSpecies;
    vect_str_t  mDonorDatabasePaths;
    vect_str_t  mRecipientDatabasePaths;
    std::string mGFFPath;
    std::string mBlastType; // string to signify blast type
    const std::string BLASTX_STR           = "blastx";
    const std::string BLASTP_STR           = "blastp";
    const fp64 DMND_TARGET_COVERAGE = 50.0;
    const fp64 DMND_QUERY_COVERAGE = 50.0;
    const fp64 DMND_E_VAL = 1E-5;

    std::string get_database_shortname(std::string &full_path);
    std::string get_database_output_path(std::string &database_name);
};


#endif //ENTAP_ABSTRACTHORIZONTALGENETRANSFER_H
