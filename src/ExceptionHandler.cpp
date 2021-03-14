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


//*********************** Includes *****************************

#include "ExceptionHandler.h"
//**************************************************************

/**
 * ======================================================================
 * Function ExceptionHandler::ExceptionHandler(const std::string& msg, int err)
 *
 * Description          - Constructor used for ExceptionHandler
 *                      - Sets EnTAP error code and error message
 *
 * Notes                - Constructor
 *
 * @param msg           - Error message
 * @param err           - Error code
 *
 * @return              - ExceptionHandler object
 *
 * =====================================================================
 */
ExceptionHandler::ExceptionHandler(const std::string& msg, int err) {
    this->mErrCode = err;
    this->mMessage = msg;
}


/**
 * ======================================================================
 * Function ExceptionHandler::print_msg()
 *
 * Description          - Prints custom error message for each error code
 *                        as well as set error message from constructor
 *
 * Notes                - None
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void ExceptionHandler::print_msg(FileSystem* filesystem) {

    std::stringstream added_msg;
    std::string out_msg;

    added_msg << "Error code: " << mErrCode << "\n";

    switch (mErrCode) {
        case ERR_ENTAP_CONFIG_CREATE:
            added_msg << "Error in creating the EnTAP configuration file. If this persists, download "
                    "the file from GitHub";
            break;
        case ERR_ENTAP_CONFIG_CREATE_SUCCESS:
            added_msg << "Configuration file was not found and was generated for you, make sure to "
                    "check the paths before continuing.";
            break;
        case ERR_ENTAP_INIT_TAX_DOWN:
            added_msg << "Error in downloading the taxonomic database."; // No longer used
            break;
        case ERR_ENTAP_INIT_TAX_INDEX:
            added_msg << "Error parsing the downloaded taxonomic database. Ensure that it was downloaded "
                      "correctly. ";
            break;
        case ERR_ENTAP_INIT_TAX_SERIAL:
            added_msg << "Error in indexing the taxonomic database. This process requires Boost Serialization "
                    "libraries as well as a properly downloaded taxonomic datbase.";
            break;
        case ERR_ENTAP_INIT_EGGNOG:
            added_msg << "Error in downloading EggNOG databases through EggNOG Python script. Ensure that you "
                    "have a proper EggNOG mapper and Python installation.";
            break;
        case ERR_ENTAP_RUN_EGGNOG:
            added_msg << "Error in running EggNOG Emapper. EggNOG requires a sqlite module in your"
                    "distribution of Python as well as a global DIAMOND installation to call from.";
            break;
        case ERR_ENTAP_RUN_EGGNOG_DMND:
            break;
        case ERR_ENTAP_PARSE_EGGNOG_DMND:
            added_msg << "Error in parsing data associated with EggNOG results.";
            break;
        case ERR_ENTAP_PARSE_EGGNOG:
            added_msg << "Error in parsing EggNOG data. Ensure that EggNOG ran properly and the output "
                      "has the proper data contained within.";
            break;
        case ERR_ENTAP_INIT_TAX_READ:
            added_msg << "Unable to read the Taxonomic database into memory, ensure the paths are "
                    "correct along with how you configured the file. You may need to remove it and "
                    "re-download.";
            break;
        case ERR_ENTAP_INIT_GO_DOWNLOAD:
            added_msg << "Error in downloading the Gene Ontology database, ensure you are connected "
                    "to the internet as well as having the Unix module wget available (this is "
                    "available by default on most systems).";
            break;
        case ERR_ENTAP_INIT_GO_UNZIP:
            added_msg << "Error in unzipping the Gene Ontology data downloaded. This could be due to "
                    "not having gzip available on your system (generally available by default). You may "
                    "want to try using the pre-configured database in the Git repo.";
            break;
        case ERR_ENTAP_INIT_GO_INDEX:
            added_msg << "Error in indexing Gene Ontology database. This could be due to a poor"
                    " download of the data or an issue with the required Boost library (serialization).";
            break;
        case ERR_ENTAP_RUN_GENEMARK_PARSE:
            added_msg << "Ensure that your GeneMarkS-T run completed successfully and you have "
                    "sequences in your output file!";
            break;
        case ERR_ENTAP_RUN_GENEMARK_STATS:
            added_msg << "Ensure the GeneMarkS-T run finished successfully. If it has, re-run"
                    " with the --overwrite flag or delete your frame_selection directory.";
            break;
        case ERR_ENTAP_RUN_GENEMARK_MOVE:
            added_msg << "Ensure GeneMarkS-T ran properly and the files are all located within "
                    "the frame_selection directory.";
            break;
        case ERR_ENTAP_RUN_RSEM_VALIDATE:
            added_msg << "Ensure that you have specified the correct path to the RSEM "
                    "directory and have it properly compiled.";
            break;
        case ERR_ENTAP_RUN_RSEM_EXPRESSION:
            added_msg << "Ensure that you have specified the correct path to the RSEM "
                    "directory and have it properly compiled.";
            break;
        case ERR_ENTAP_RUN_RSEM_EXPRESSION_PARSE:
            added_msg << "Ensure that RSEM ran properly and the output has sequences in it!";
            break;
        case ERR_ENTAP_RUN_SIM_SEARCH_FILTER:
            added_msg << "Ensure the similarity searching finished properly and your output files"
                    " are not empty.";
            break;
        case ERR_ENTAP_RUN_SIM_SEARCH_RUN:
            added_msg << "Ensure you have specified the proper path to the executable. Additionally, "
                    "the database you are specifying is configured for DIAMOND (.dmnd).";
            break;
        case ERR_ENTAP_RUN_ANNOTATION:
            added_msg << "Ensure you have specified the proper path to the executable.";
            break;
        case ERR_ENTAP_DATABASE_QUERY:
            added_msg << "Ensure that your EggNOG run finished successfully and the files are not empty.";
            break;
        case ERR_ENTAP_RUN_INTERPRO:
            added_msg << "Ensure you have specified the proper path to the InterProScan executable "
                    "and it is properly compiled on your system. Additionally, ensure you have the "
                    "proper databases downloaded that you would like to run against.";
            break;
        case ERR_ENTAP_INIT_GEN_SERIAL_DATA:
            added_msg << "Ensure you have the proper Boost libraries installed as well as an " \
                         "internet connection. You may also try downloading them!";
            break;
        case ERR_ENTAP_INIT_GEN_SQL_DATA:
            added_msg << "Ensure you have a proper internet connection";
            break;
        case ERR_ENTAP_INIT_DOWN_SERIAL_DATA:
            added_msg << "Ensure you have a proper internet connection and wget available";
            break;
        case ERR_ENTAP_INIT_DOWN_SQL_DATA:
            added_msg << "Ensure you have a proper internet connection and wget available";
            break;
        default:
            FS_dprint("Error code (" + std::to_string(mErrCode) + ") not recognized");
//            added_msg << "Error code not recognized.";
            break;
    }

    added_msg << "\n" << what();
    out_msg = added_msg.str();
    if (filesystem != nullptr) FS_dprint(out_msg);
    std::cerr << out_msg << std::endl;
}


const char* ExceptionHandler::what() {
    return mMessage.c_str();
}


uint16 ExceptionHandler::getErr_code() const {
    return mErrCode;
}