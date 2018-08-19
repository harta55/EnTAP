/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
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

#ifndef ENTAP_MODEGGNOGDMND_H
#define ENTAP_MODEGGNOGDMND_H


#include "AbstractOntology.h"

class ModEggnogDMND : public AbstractOntology {

public:
    ModEggnogDMND(std::string &out, std::string &in_hits, std::string &ont_out, bool blastp, std::vector<uint16> &lvls,
                  EntapDataPtrs &entap_data, std::string sql_db_path);
    virtual std::pair<bool, std::string> verify_files() override;
    ~ModEggnogDMND();
    virtual void execute() override ;
    virtual void parse() override;
    static bool is_executable();

private:
    std::string get_output_dmnd_filepath(bool final);

    static constexpr int DMND_COL_NUMBER = 14;
    const std::string EGGNOG_DMND_DIR = "EggNOG_DMND/";
    std::string _out_hits;
    std::string _egg_out_dir;
    std::string _eggnog_db_path;


};


#endif //ENTAP_MODEGGNOGDMND_H
