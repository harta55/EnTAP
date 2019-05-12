/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
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

AbstractOntology::AbstractOntology(std::string &in_hits, std::string &ont_out, EntapDataPtrs &entap_data,
                                   std::string mod_name, std::string &exe)
: EntapModule(ont_out, in_hits, entap_data, mod_name, exe) {

    _go_levels          = _pUserInput->get_user_input<vect_uint16_t>(_pUserInput->INPUT_FLAG_GO_LEVELS);
    _execution_state    = GENE_ONTOLOGY;
}