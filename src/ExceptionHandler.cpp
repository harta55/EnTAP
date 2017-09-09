/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/


//*********************** Includes *****************************
#include <iostream>
#include <ctime>
#include <fstream>
#include <exception>
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
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
    this->err_code = err;
    this->message = msg;
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
void ExceptionHandler::print_msg() {

    std::stringstream added_msg;
    std::string out_msg;

    added_msg << "Error code: " << err_code << "\n";

    switch (err_code) {
        case ENTAP_ERR::E_INPUT_PARSE:
            added_msg << "Error in parsing input data, please consult -h for more information.";
            break;
        case ENTAP_ERR::E_CONFIG_PARSE:
            added_msg << "Error in parsing the EnTAP configuration file, ensure all parameters are"
                    "in the correct format.";
            break;
        case ENTAP_ERR::E_CONFIG_CREATE:
            added_msg << "Error in creating the EnTAP configuration file. If this persists, download"
                    "the file from GitHub";
            break;
        case ENTAP_ERR::E_CONFIG_CREATE_SUCCESS:
            added_msg << "Configuration file was not found and was generated, make sure to "
                    "check the paths before continuing.";
            break;
        case ENTAP_ERR::E_INIT_TAX_DOWN:
            added_msg << "Error in downloading the taxonomic database. Ensure that the script is at "
                      << TAX_DOWNLOAD_EXE;
            break;
        case ENTAP_ERR::E_INIT_TAX_INDEX:
            added_msg << "Error parsing the downloaded taxonomic database. Ensure that it was downloaded "
                      "correctly. ";
            break;
        case ENTAP_ERR::E_INIT_TAX_SERIAL:
            added_msg << "Error in indexing the taxonomic database. This process requires Boost Serialization "
                    "libraries as well as a properly downloaded taxonomic datbase.";
            break;
        case ENTAP_ERR::E_INIT_EGGNOG:
            added_msg << "Error in downloading EggNOG databases through EggNOG Python script. Ensure that you "
                    "have a proper EggNOG mapper and Python installation.";
            break;
        case ENTAP_ERR::E_RUN_EGGNOG:
            added_msg << "Error in running EggNOG Emapper. EggNOG requires a sqlite module in your"
                    "distribution of Python as well as a global DIAMOND installation to call from.";
            break;
        case ENTAP_ERR::E_PARSE_EGGNOG:
            added_msg << "Error in parsing EggNOG data. Ensure that EggNOG ran properly and the output "
                      "has the proper data contained within.";
            break;
        default:
            added_msg << "Error code not recognized.";
            break;
    }

    added_msg << "\n" << what();
    out_msg = added_msg.str();
    print_debug(out_msg);
    std::cerr << out_msg << std::endl;
}


const char* ExceptionHandler::what() {
    return message.c_str();
}


int ExceptionHandler::getErr_code() const {
    return err_code;
}




