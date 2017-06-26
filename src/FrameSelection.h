//
// Created by harta on 5/7/17.
//

#ifndef ENTAP_FRAMESELECTION_H
#define ENTAP_FRAMESELECTION_H
#include <iostream>
#include <map>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"


class FrameSelection {
    struct frame_seq {
        unsigned long length;
        std::string sequence;
        std::string frame_type;
    };
    typedef std::map<std::string,FrameSelection::frame_seq> frame_map_type;

    public:
        std::string execute(std::string,std::map<std::string,QuerySequence>&);
        FrameSelection(std::string&, std::string&, std::string&,
                       boost::program_options::variables_map &);

    private:
        std::string _exe_path,_inpath, _outpath;
        bool _overwrite;
        short _software_flag;
        std::string genemarkst(std::map<std::string,QuerySequence> &);
        void genemarkStats(std::string&,std::string&,std::map<std::string,QuerySequence> &);
        frame_map_type genemark_parse_protein(std::string&);
        void genemark_parse_lst(std::string &, std::map<std::string, frame_seq>&);
};


#endif //ENTAP_FRAMESELECTION_H
