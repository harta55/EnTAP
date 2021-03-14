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


#ifndef ENTAP_ABSTRACTONTOLOGY_H
#define ENTAP_ABSTRACTONTOLOGY_H

#include "../common.h"
#include "../EntapGlobals.h"
#include "../EntapModule.h"

class AbstractOntology : public EntapModule {
public:
    AbstractOntology(std::string &in_hits, std::string &ont_out,
                     EntapDataPtrs &entap_data, std::string mod_name,
                     std::vector<ENTAP_HEADERS> &module_headers);
    ~AbstractOntology() = default;
    virtual ModVerifyData verify_files()=0;
    virtual void execute() = 0;
    virtual void parse() = 0;
    virtual bool set_version() = 0;

    virtual void set_success_flags() override ;
};


#endif //ENTAP_ABSTRACTONTOLOGY_H
