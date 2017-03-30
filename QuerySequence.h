//
// Created by harta on 3/29/17.
//

#ifndef ENTAP_QUERYSEQUENCE_H
#define ENTAP_QUERYSEQUENCE_H
#include <iostream>


class QuerySequence {
    public:
        QuerySequence operator==(const QuerySequence& querySequence);
        QuerySequence(float);
        QuerySequence();

    private:
        float e_val;
        std::string database_path;



};


#endif //ENTAP_QUERYSEQUENCE_H
