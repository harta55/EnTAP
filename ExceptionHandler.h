//
// Created by Alex on 9/14/2016.
//

#ifndef ENTAP_EXCEPTIONHANDLER_H
#define ENTAP_EXCEPTIONHANDLER_H

#include <string>

class ExceptionHandler {
    public:
        ExceptionHandler(const int err_tag, std::string m);
        int getErrTag();
        static const short int except_help = 0;
        static const short int except_input_parse = 1;

    private:
        int err_tag;
        void parse_flag(int flag, std::string msg);
        void print_msg(std::string, bool b);
    };

#endif //ENTAP_EXCEPTIONHANDLER_H
