//
// Created by Alex on 9/14/2016.
//

#ifndef ENTAP_EXCEPTIONHANDLER_H
#define ENTAP_EXCEPTIONHANDLER_H

#include <string>
#include <exception>
#include <sstream>

class ExceptionHandler: public std::runtime_error{
    int err;
    std::string message;
    public:
        ExceptionHandler(const std::string &__arg, const int err_tag, std::string m);
        virtual const char *what() const noexcept override ;

    int getErrTag();
        static const short int except_help = 0;
        static const short int except_input_parse = 1;

    private:
        static std::ostringstream cnvt;
        int err_tag;
        void parse_flag(int flag, std::string msg);
        void print_msg(std::string, bool b);
    };

#endif //ENTAP_EXCEPTIONHANDLER_H
