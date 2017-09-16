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


#ifndef ENTAP_MODEGGNOG_H
#define ENTAP_MODEGGNOG_H


#include "AbstractOntology.h"

class ModEggnog : public AbstractOntology{

public:
    ModEggnog(std::string &exe,std::string &out, std::string &in,
              std::string &in_no_hits,std::string &proc,
            std::string &fig, std::string &ont, GraphingManager *graphing) :
    AbstractOntology(exe, out, in, in_no_hits,proc, fig, ont, graphing){}

    virtual std::pair<bool, std::string> verify_files() override ;
    virtual void execute(std::map<std::string, QuerySequence>&) override ;
    virtual void parse(std::map<std::string, QuerySequence>&) override ;
    virtual void set_data(std::vector<short>&, std::string&, int) override ;

private:

    static constexpr short EGGNOG_COL_NUM     = 12;
    const std::string GRAPH_EGG_TAX_BAR_TITLE = "Top_10_Tax_Levels";
    const std::string GRAPH_EGG_TAX_BAR_PNG   = "eggnog_tax_scope.png";
    const std::string GRAPH_EGG_TAX_BAR_TXT   = "eggnog_tax_scope.txt";

    int _threads;
    std::string _eggnog_db_path;
    std::string _out_no_hits;
    std::string _out_hits;
    std::vector<short> _go_levels;

    std::string eggnog_format(std::string);
};


#endif //ENTAP_MODEGGNOG_H
