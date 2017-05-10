//
// Created by harta on 5/10/17.
//

#ifndef ENTAP_SIMILARITYSEARCH_H
#define ENTAP_SIMILARITYSEARCH_H
#include <iostream>
#include <list>


class SimilaritySearch {
public:
    std::list<std::string> execute(short, std::string, bool);
    SimilaritySearch(std::list<std::string>&, std::string, int, bool, std::string,std::string);
private:
    std::list<std::string> _database_paths;
    std::string _exe, _outpath, _input;
    int _threads; bool _overwrite,_blastp;
    std::list<std::string> diamond();
    void diamond_blast(std::string, std::string, std::string,std::string&,int&, std::string&);
    std::list<std::string> verify_diamond_files(std::string&,
        std::string&, std::string);
};


#endif //ENTAP_SIMILARITYSEARCH_H
