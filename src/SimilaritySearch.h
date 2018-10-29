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

#ifndef ENTAP_SIMILARITYSEARCH_H
#define ENTAP_SIMILARITYSEARCH_H

//*********************** Includes *****************************
#include "common.h"
#include "QuerySequence.h"
#include "GraphingManager.h"
#include "QueryData.h"
#include "database/EntapDatabase.h"
#include "FileSystem.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
#include "UserInput.h"
#include "database/EntapDatabase.h"

//**************************************************************


class SimilaritySearch {

public:

    //******************** Public Prototype Functions *********************
    std::vector<std::string> execute(std::string, bool);
    SimilaritySearch(databases_t&, std::string,EntapDataPtrs&);
    SimilaritySearch();
    void parse_files(std::string);
    static bool is_executable();
    //**************************************************************


private:

    const uint8       UNIPROT_ATTEMPTS                           = 15;   // Number of attempts to see if database is uniprot
    const std::string _NCBI_REGEX                                = "\\[(.+)\\](?!.+\\[.+\\])";
    const std::string _UNIPROT_REGEX                             = "OS=(.+?)\\s\\S\\S=";
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
    const std::string SIM_SEARCH_DIR                             = "similarity_search/";
    const std::string PROCESSED_DIR                              = "processed/";
    const std::string RESULTS_DIR                                = "overall_results/";
    const std::string FIGURE_DIR                                 = "figures/";

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

    const std::string BLASTX                                     = "blastx";
    const std::string BLASTP                                     = "blastp";


    static constexpr int DMND_COL_NUMBER = 14;
    static constexpr short COUNT_TOP_SPECIES = 20;

    const std::vector<const std::string*> DEFAULT_HEADERS {
            &ENTAP_EXECUTE::HEADER_QUERY,
            &ENTAP_EXECUTE::HEADER_SUBJECT,
            &ENTAP_EXECUTE::HEADER_PERCENT,
            &ENTAP_EXECUTE::HEADER_ALIGN_LEN,
            &ENTAP_EXECUTE::HEADER_MISMATCH,
            &ENTAP_EXECUTE::HEADER_GAP_OPEN,
            &ENTAP_EXECUTE::HEADER_QUERY_S,
            &ENTAP_EXECUTE::HEADER_QUERY_E,
            &ENTAP_EXECUTE::HEADER_SUBJ_S,
            &ENTAP_EXECUTE::HEADER_SUBJ_E,
            &ENTAP_EXECUTE::HEADER_E_VAL,
            &ENTAP_EXECUTE::HEADER_COVERAGE,
            &ENTAP_EXECUTE::HEADER_TITLE,
            &ENTAP_EXECUTE::HEADER_SPECIES,
            &ENTAP_EXECUTE::HEADER_DATABASE,
            &ENTAP_EXECUTE::HEADER_FRAME,
            &ENTAP_EXECUTE::HEADER_CONTAM,
            &ENTAP_EXECUTE::HEADER_INFORM,
            &ENTAP_EXECUTE::HEADER_UNI_DATA_XREF,
            &ENTAP_EXECUTE::HEADER_UNI_COMMENTS
    };

    std::vector<std::string>        _database_paths;
    std::vector<std::string>        _sim_search_paths;
    std::vector<std::string>        _uninformative_vect;
    std::string                     _diamond_exe;
    std::string                     _outpath;
    std::string                     _input_path;
    std::string                     _sim_search_dir;
    std::string                     _processed_path;
    std::string                     _figure_path;
    std::string                     _results_path;
    std::string                     _input_lineage;
    std::string                     _input_species;
    std::string                     _blast_type;
    std::string                     _transcript_shortname;
    int                             _threads;
    bool                            _overwrite;
    bool                            _blastp;
    fp64                            _e_val;
    fp32                            _qcoverage;
    fp32                            _tcoverage;
    uint8                           _software_flag;
    std::vector<std::string>        _contaminants;
    GraphingManager                 *_pGraphingManager;
    QueryData                       *_pQUERY_DATA;
    FileSystem                      *_pFileSystem;
    UserInput                       *_pUserInput;
    EntapDatabase                   *_pEntapDatabase;
    std::unordered_map<std::string,std::string> _file_to_database;

    std::vector<std::string> diamond();
    void diamond_blast(std::string, std::string, std::string,std::string&,int&, std::string&);
    std::vector<std::string> verify_diamond_files();
    void diamond_parse(std::vector<std::string>&);
    std::pair<bool,std::string> is_contaminant(std::string, std::vector<std::string>&);
    bool is_informative(std::string);
    bool is_uniprot_entry(std::string& accession, UniprotEntry &entry);
    void print_header(std::ofstream&);
    std::string get_species(std::string &title);
    void calculate_best_stats (bool,std::string="");
    std::string get_database_shortname(std::string&);
    std::string get_transcriptome_shortname();

};


#endif //ENTAP_SIMILARITYSEARCH_H
