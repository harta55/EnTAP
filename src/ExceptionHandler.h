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

#ifndef ENTAP_EXCEPTIONHANDLER_H
#define ENTAP_EXCEPTIONHANDLER_H

#include "EntapGlobals.h"
#include "FileSystem.h"

// TODO fix/update error code handling
enum ENTAP_ERR {

    ERR_ENTAP_SUCCESS                     = 0u,
    ERR_ENTAP_INPUT_PARSE                 = 10u,
    ERR_ENTAP_CONFIG_PARSE                = 12u,
    ERR_ENTAP_CONFIG_CREATE               = 13u,
    ERR_ENTAP_CONFIG_CREATE_SUCCESS       = 14u,
    ERR_ENTAP_INIT_TAX_DOWN               = 20u,
    ERR_ENTAP_INIT_TAX_INDEX              = 21u,
    ERR_ENTAP_INIT_TAX_SERIAL             = 22u,
    ERR_ENTAP_INIT_DOWNLOAD               = 23u,
    ERR_ENTAP_INIT_INDX_DATA_NOT_FOUND    = 30u,
    ERR_ENTAP_INIT_INDX_DATABASE          = 31u,
    ERR_ENTAP_INIT_GEN_SERIAL_DATA        = 32u,
    ERR_ENTAP_INIT_GEN_SQL_DATA           = 33u,
    ERR_ENTAP_INIT_DOWN_SQL_DATA          = 34u,
    ERR_ENTAP_INIT_DOWN_SERIAL_DATA       = 35u,
    ERR_ENTAP_INIT_DATA_GENERIC           = 36u,
    ERR_ENTAP_INIT_BUSCO_GENERIC          = 37u,
    ERR_ENTAP_INIT_EGGNOG                 = 40u,

    ERR_ENTAP_INIT_TAX_READ               = 55u,
    ERR_ENTAP_INIT_GO_DOWNLOAD            = 60u,
    ERR_ENTAP_INIT_GO_UNZIP               = 61u,
    ERR_ENTAP_INIT_GO_PARSE               = 62u,
    ERR_ENTAP_INIT_GO_INDEX               = 63u,
    ERR_ENTAP_READ_ENTAP_SERIAL_DATA      = 70u,
    ERR_ENTAP_READ_ENTAP_SQL_DATA         = 71u,
    ERR_ENTAP_READ_ENTAP_DATA_GENERIC     = 72u,
    ERR_ENTAP_GENERATE_TRANSCRIPTOME      = 80u,

    ERR_ENTAP_RUN_GENEMARK                = 100u,
    ERR_ENTAP_RUN_GENEMARK_PARSE          = 101u,
    ERR_ENTAP_RUN_GENEMARK_STATS          = 102u,
    ERR_ENTAP_RUN_GENEMARK_MOVE           = 103u,
    ERR_ENTAP_RUN_TRANSDECODER            = 105u,
    ERR_ENTAP_RUN_TRANSDECODER_PARSE      = 106u,
    ERR_ENTAP_RUN_TRANSDECODER_STATS      = 107u,
    ERR_ENTAP_RUN_TRANSDECODER_MOVE       = 108u,
    ERR_ENTAP_RUN_RSEM_VALIDATE           = 110u,
    ERR_ENTAP_RUN_RSEM_CONVERT            = 111u,
    ERR_ENTAP_RUN_RSEM_EXPRESSION         = 112u,
    ERR_ENTAP_RUN_RSEM_EXPRESSION_PARSE   = 113u,
    ERR_ENTAP_RUN_FILTER                  = 120u,
    ERR_ENTAP_RUN_SIM_SEARCH_FILTER       = 140u,
    ERR_ENTAP_RUN_SIM_SEARCH_RUN          = 141u,
    ERR_ENTAP_RUN_ANNOTATION              = 150u,
    ERR_ENTAP_GO_READ                     = 151u,
    ERR_ENTAP_RUN_EGGNOG                  = 160u,
    ERR_ENTAP_DATABASE_QUERY              = 161u,
    ERR_ENTAP_PARSE_EGGNOG                = 162u,
    ERR_ENTAP_EGGNOG_FILES                = 163u,   // Some EggNOG file missing
    ERR_ENTAP_RUN_INTERPRO                = 170u,
    ERR_ENTAP_PARSE_INTERPRO              = 171u,
    ERR_ENTAP_PARSE_EGGNOG_DMND           = 180u,
    ERR_ENTAP_RUN_EGGNOG_DMND             = 181u,
    ERR_ENTAP_RUN_BUSCO                   = 190u,
    ERR_ENTAP_PARSE_BUSCO                 = 191u,
    ERR_ENTAP_VERSION_UNSUPPORTED_BUSCO   = 192u,
    ERR_ENTAP_FILE_IO                     = 200u,
    ERR_ENTAP_MEM_ALLOC                   = 201u,
    ERR_ENTAP_MAX                         = 201u
};


class ExceptionHandler: public std::exception{

public:
    ExceptionHandler(const std::string&, int);
    const char* what();
    uint16 getErr_code() const;
    void print_msg(FileSystem*);

private:
    uint16 mErrCode;
    std::string mMessage;

//    static const std::string ERR_ENTAP_STR [ERR_ENTAP_MAX];
};

#endif //ENTAP_EXCEPTIONHANDLER_H
