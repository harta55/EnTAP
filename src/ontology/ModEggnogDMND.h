/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
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
#include "../config.h"

class ModEggnogDMND : public AbstractOntology {

public:
    ModEggnogDMND(std::string &ont_out, std::string &in_hits,
                  EntapDataPtrs &entap_data, std::string &exe, std::string sql_db_path);
    virtual ModVerifyData verify_files() override;
    ~ModEggnogDMND();
    virtual void execute() override ;
    virtual void parse() override;
    static bool is_executable(std::string &exe);

    static const std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;

private:
    std::string get_output_dmnd_filepath(bool final);
    void calculate_stats(std::stringstream &stream);

    static constexpr int DMND_COL_NUMBER = 14;
    const uint32      STATUS_UPDATE_HITS = 5000;
    const std::string GRAPH_EGG_TAX_BAR_TITLE = "Top_Tax_Levels";
    const std::string GRAPH_EGG_TAX_BAR_PNG   = "eggnog_tax_scope.png";
    const std::string GRAPH_EGG_TAX_BAR_TXT   = "eggnog_tax_scope.txt";
    const std::string EGG_ANNOT_RESULTS       = "annotation_results";
    const std::string EGG_ANNOT_STD           = "annotation_std";
    const std::string EGG_ANNOT_APPEND        = ".emapper.annotations";
    static constexpr uint16 COUNT_TOP_TAX_SCOPE = 10;

    std::string _out_hits;
    std::string _eggnog_db_path;
};


#endif //ENTAP_MODEGGNOGDMND_H
