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

#include "AbstractExpression.h"
#include "../QueryData.h"

AbstractExpression::AbstractExpression(std::string &execution_stage_path, std::string &in_hits,
                                       EntapDataPtrs &entap_data, std::string module_name,
                                       std::vector<ENTAP_HEADERS> &module_headers) :
EntapModule(execution_stage_path, in_hits, entap_data, module_name, module_headers) {

    if (mpUserInput->has_input(INPUT_FLAG_ALIGN)) { // Will be true
        mAlignPath = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_ALIGN);
    }
    mExecutionState= EXPRESSION_FILTERING;
    mIsSingle      = mpUserInput->has_input(INPUT_FLAG_SINGLE_END);
    mFPKM          = mpUserInput->get_user_input<ent_input_fp_t >(INPUT_FLAG_FPKM);
}

void AbstractExpression::set_success_flags() {
    mpQueryData->set_is_success_expression(true);
}
