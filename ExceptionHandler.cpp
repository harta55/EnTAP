#include <iostream>
#include <ctime>
#include <fstream>
#include <exception>
#include "ExceptionHandler.h"


static std::ostringstream cnvt;

ExceptionHandler::ExceptionHandler(const std::string &__arg, const int err_tag, std::string m) : runtime_error(__arg) {
    err = err_tag;
    message = m;
}

int ExceptionHandler::getErrTag() {
    return err_tag;
}

void ExceptionHandler::parse_flag(int flag, std::string msg) {
    switch (flag) {
        case 0:
            // help flag, nothing to print
            return;
        case 1:
            // error in input
            print_msg("Error in parsing input data, please check with -h (help) "
                              "and try again.", true);
            if (msg.compare("") != 0){
                print_msg("Error: "+ msg, true);
            }
            return;
        default:
            return;
    }
}

void ExceptionHandler::print_msg(std::string msg, bool b) {
    time_t rawtime;
    time(&rawtime);
    std::string date_time = ctime(&rawtime);
    std::ofstream log_file(
            "debug.txt", std::ios_base::out | std::ios_base::app );
    if (b) {
        //error
        log_file << date_time.substr(0, date_time.size() - 2)
                     + ": " + msg << std::endl;
    }else {
        log_file << date_time.substr(0, date_time.size() - 2)
                     + ": " + msg << std::endl;
    }
}

virtual const char* what() const throw()
{
    cnvt.str( "" );

    cnvt << std::runtime_error::what();
    return cnvt.str().c_str();
}




