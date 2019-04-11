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

#ifndef ENTAP_MODDIAMOND_H
#define ENTAP_MODDIAMOND_H


#include "AbstractSimilaritySearch.h"

class ModDiamond : public AbstractSimilaritySearch {

public:
    ModDiamond(std::string &out, std::string &in_hits,EntapDataPtrs &entap_data,
                std::string &exe, vect_str_t &databases);
    ~ModDiamond() override = default;

    // ModEntap overrides
    virtual ModVerifyData verify_files() override;
    virtual void execute() override ;
    virtual void parse() override ;
    static bool is_executable(std::string& exe);

    // AbstractSimilaritySearch overrides
    virtual bool run_blast(SimSearchCmd *cmd, bool use_defaults);

    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;


private:
    static constexpr int DMND_COL_NUMBER = 14;
    const std::string SIM_SEARCH_DATABASE_BEST_TSV               = "best_hits.tsv";
    const std::string SIM_SEARCH_DATABASE_BEST_TSV_NO_CONTAM     = "best_hits_no_contam.tsv";
    const std::string SIM_SEARCH_DATABASE_BEST_FA_NUCL           = "best_hits.fnn";
    const std::string SIM_SEARCH_DATABASE_BEST_FA_PROT           = "best_hits.faa";
    const std::string SIM_SEARCH_DATABASE_BEST_FA_NUCL_NO_CONTAM = "best_hits_no_contam.fnn";
    const std::string SIM_SEARCH_DATABASE_BEST_FA_PROT_NO_CONTAM = "best_hits_no_contam.faa";

    const std::string SIM_SEARCH_DATABASE_CONTAM_TSV             = "best_hits_contam.tsv";
    const std::string SIM_SEARCH_DATABASE_CONTAM_FA_NUCL         = "best_hits_contam.fnn";
    const std::string SIM_SEARCH_DATABASE_CONTAM_FA_PROT         = "best_hits_contam.faa";
    const std::string SIM_SEARCH_DATABASE_NO_HITS_NUCL           = "no_hits.fnn";
    const std::string SIM_SEARCH_DATABASE_NO_HITS_PROT           = "no_hits.faa";
    const std::string SIM_SEARCH_DATABASE_UNSELECTED             = "unselected.tsv";

    // Graphing constants
    const uint8 GRAPH_SOFTWARE_FLAG                              = 3;
    const uint8 GRAPH_BAR_FLAG                                   = 1;
    const uint8 GRAPH_SUM_FLAG                                   = 2;
    const std::string GRAPH_DATABASE_SUM_TITLE                   = "_Summary";
    const std::string GRAPH_DATABASE_SUM_TXT                     = "_summary_bar.txt";
    const std::string GRAPH_DATABASE_SUM_PNG                     = "_summary_bar.png";
    const std::string GRAPH_SPECIES_BAR_TXT                      = "_species_bar.txt";
    const std::string GRAPH_SPECIES_BAR_PNG                      = "_species_bar.png";
    const std::string GRAPH_SPECIES_TITLE                        = "_Top_10_Species_Distribution";
    const std::string GRAPH_CONTAM_BAR_TXT                       = "_contam_bar.txt";
    const std::string GRAPH_CONTAM_BAR_PNG                       = "_contam_bar.png";
    const std::string GRAPH_CONTAM_TITLE                         = "_Top_10_Contaminant_Distribution";
    const std::string UNINFORMATIVE_FLAG                         = "Uninformative";
    const std::string INFORMATIVE_FLAG                           = "Informative";
    const std::string NO_HIT_FLAG                                = "No Hits";

    void calculate_best_stats(bool is_final, std::string database_path="");
};


#endif //ENTAP_MODDIAMOND_H
