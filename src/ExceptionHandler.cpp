#include <iostream>
#include <ctime>
#include <fstream>
#include <exception>
#include "ExceptionHandler.h"
#include "EntapGlobals.h"


ExceptionHandler::ExceptionHandler(const std::string& msg, int err) {
    this->err_code = err;
    this->message = msg;
}

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




