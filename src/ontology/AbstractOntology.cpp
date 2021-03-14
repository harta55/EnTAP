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


#include "AbstractOntology.h"
#include "../QueryData.h"

AbstractOntology::AbstractOntology(std::string &in_hits, std::string &ont_out, EntapDataPtrs &entap_data,
                                   std::string mod_name, std::vector<ENTAP_HEADERS> &module_headers)
: EntapModule(ont_out, in_hits, entap_data, mod_name, module_headers) {

    mExecutionState    = GENE_ONTOLOGY;
}

void AbstractOntology::set_success_flags() {
    mpQueryData->set_is_success_ontology(true);
}
