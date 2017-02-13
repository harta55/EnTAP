#include <iostream>
#include <ctime>
#include <fstream>
#include <exception>
#include "ExceptionHandler.h"
#include "ErrorFlags.h"


ExceptionHandler::ExceptionHandler(const std::string& msg, int err) {
    this->err_code = err;
    this->message = msg;
}

//void ExceptionHandler::parse_flag(int flag, std::string msg) {
//    switch (flag) {
//        case 0:
//            // help flag, nothing to print
//            return;
//        case 1:
//            // error in input
//            print_msg("Error in parsing input data, please check with -h (help) "
//                              "and try again.", true);
//            if (msg.compare("") != 0){
//                print_msg("Error: "+ msg, true);
//            }
//            return;
//        default:
//            return;
//    }
//}

void ExceptionHandler::print_msg() {
    time_t rawtime;
    time(&rawtime);
    std::string date_time = ctime(&rawtime);
    std::ofstream log_file(
            "debug.txt", std::ios_base::out | std::ios_base::app );
    std::string added_msg;

    if (err_code == ENTAPERR::E_INPUT_PARSE) std::cout<<"ERROR!"<<std::endl;
    switch (err_code) {
        case ENTAPERR::E_INPUT_PARSE:
            added_msg = "Error in parsing input data, please consult -h for more "
                    "information.";
            break;
        default:
            added_msg = "Error code not recognized.";
    }
    log_file << date_time.substr(0, date_time.size() - 2)
                 + ": " + added_msg << std::endl ;
    log_file <<what()<<std::endl;
}

const char* ExceptionHandler::what() {
    return message.c_str();
}

int ExceptionHandler::getErr_code() const {
    return err_code;
}




