/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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

#ifndef ENTAP_BASESIMILARITYSEARCH_H
#define ENTAP_BASESIMILARITYSEARCH_H


#include "../EntapModule.h"

class AbstractSimilaritySearch : public EntapModule{

public:

    struct SimSearchCmd {
        std::string         exe_path;
        std::string         query_path;
        std::string         output_path;
        std::string         database_path;
        std::string         std_out_path;

        bool                blastp;
        bool                more_sensitive;

        uint16              threads;
        fp64                eval;
        fp64                qcoverage;      // query coverage
        fp64                tcoverage;      // target coverage
        uint16              top_num;        // default = 3
        uint16              output_flags;    // currently unused
    };

    AbstractSimilaritySearch(std::string &execute_stage_path,
                             std::string &in_hits,
                             EntapDataPtrs &entap_data,
                             std::string mod_name,
                             std::vector<ENTAP_HEADERS> &module_headers,
                             vect_str_t &databases);
    ~AbstractSimilaritySearch() = default;
    virtual ModVerifyData verify_files()=0;
    virtual void execute() = 0;
    virtual void parse() = 0;
    virtual void set_success_flags() override ;
    virtual bool set_version() = 0;

    virtual bool run_blast(SimSearchCmd *cmd, bool use_defaults) = 0;

protected:

    vect_str_t                      mDatabasePaths;
    vect_str_t                      mUninformativeTags;
    vect_str_t                      mContaminateTaxons;
    vect_str_t                      mOutputPaths;
    std::map<std::string, std::string> mPathToDatabase;      // mapping of full database file path to shortened name
    std::string                     mInputLineage;
    std::string                     mInputSpecies;
    std::string                     mBlastType;            // string to signify blast type
    fp64                            mEVal;
    fp64                            mQCoverage;
    fp64                            mTCoverage;

    const std::string SIM_SEARCH_DATABASE_BEST_HITS              = "best_hits";
    const std::string SIM_SEARCH_DATABASE_BEST_HITS_CONTAM       = "best_hits_contam";
    const std::string SIM_SEARCH_DATABASE_BEST_HITS_NO_CONTAM    = "best_hits_no_contam";
    const std::string SIM_SEARCH_DATABASE_NO_HITS                = "no_hits";
    const std::string SIM_SEARCH_DATABASE_UNSELECTED             = "unselected";

    const std::string BLASTX_STR           = "blastx";
    const std::string BLASTP_STR           = "blastp";
    const uint8       UNIPROT_ATTEMPTS     = 15;   // Number of attempts to see if database is uniprot
    const uint8       MIN_CONTAM_COUNT     = 1;    // Minimum number of contaminants to graph
    const std::string NCBI_REGEX          = "\\[(.+)\\](?!.+\\[.+\\])";
    const std::string UNIPROT_REGEX       = "OS=(.+?)\\s\\S\\S=";

    std::string get_database_shortname(std::string &full_path);
    std::string get_database_output_path(std::string &database_name);
    std::pair<bool, std::string> is_contaminant(std::string lineage, vect_str_t &contams);
    bool is_informative(std::string title, vect_str_t &uninformative_vect);
    std::string get_species(std::string &title);
};


#endif //ENTAP_BASESIMILARITYSEARCH_H
