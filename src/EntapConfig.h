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

#ifndef ENTAPCONFIG_H
#define ENTAPCONFIG_H

//*********************** Includes *****************************
#include "UserInput.h"
#include "EntapGlobals.h"
#include "ExceptionHandler.h"
#include "FileSystem.h"
#include "database/EntapDatabase.h"
//**************************************************************

namespace entapConfig {
    //****************** Global Prototype Functions******************
    void execute_main(UserInput*, FileSystem*);
    //***************************************************************
}

#endif //ENTAP_INITHANDLER_H
