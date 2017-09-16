/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
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

//
//#ifndef ENTAP_MODINTERPRO_H
//#define ENTAP_MODINTERPRO_H
//
//
//#include "AbstractOntology.h"
//
//class ModInterpro : public AbstractOntology{
//
//    struct interpro_struct {
//        double _eval;
//        std::map<std::string,std::string> _results;
//        std::map<std::string,std::vector<std::string>> _go_map;
//    };
//
//public:
//    ModInterpro(std::string &exe, std::string &entap,std::string &out, std::string &in,
//            std::string &in_no_hits,std::string &proc,
//            std::string &fig, std::string &ont, GraphingManager *graphing) :
//    AbstractOntology(exe, entap, out, in, in_no_hits,proc, fig, ont, graphing){}
//
//    virtual std::pair<bool, std::string> verify_files() override ;
//    virtual void execute(std::map<std::string, QuerySequence>&) override ;
//    virtual void parse(std::map<std::string, QuerySequence>&) override ;
//    virtual void set_data(std::vector<short>&, std::string&, int) override ;
//
//
//private:
//    static constexpr short INTERPRO_COL_NUM = 15;
//
//};
//
//
//#endif //ENTAP_MODINTERPRO_H
