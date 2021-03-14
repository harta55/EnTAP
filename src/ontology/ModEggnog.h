#ifdef EGGNOG_MAPPER

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


#ifndef ENTAP_MODEGGNOG_H
#define ENTAP_MODEGGNOG_H


#include "AbstractOntology.h"
#include "../common.h"

/**
 * ======================================================================
 * @class ModEggnog
 *
 * Description          - This EnTAP module supports execution, parsing, and
 *                        statistical analysis of the Eggnog-mapper software
 *                        through terminal commands
 *                      - EggNOG assigns functional information to transcripts
 *                        input from the user by accessing EggNOG gene family
 *                        databases
 *                      - Parsed data is added to QueryData class
 *                      - Inherits from AbstractOntology and EntapModule classes
 *
 * Citation             - eggNOG 4.5: a hierarchical orthology framework with
 *                        improved functional annotations for eukaryotic,
 *                        prokaryotic and viral sequences. Jaime Huerta-Cepas,
 *                        Damian Szklarczyk, Kristoffer Forslund, Helen Cook,
 *                        Davide Heller, Mathias C. Walter, Thomas Rattei,
 *                        Daniel R. Mende, Shinichi Sunagawa, Michael Kuhn,
 *                        Lars Juhl Jensen, Christian von Mering, and Peer Bork.
 *                        Nucl. Acids Res. (04 January 2016) 44 (D1): D286-D293.
 *                        doi: 10.1093/nar/gkv1248
 *
 * ======================================================================
 */
class ModEggnog : public AbstractOntology{

public:
    ModEggnog(std::string &out, std::string &in, std::string &ont,
            bool blastp,std::vector<uint16>& lvls, EntapDataPtrs &entap_data, std::string&);

    ~ModEggnog();
    virtual std::pair<bool, std::string> verify_files() override ;
    virtual void execute() override ;
    virtual void parse() override ;
    static bool is_executable() ;
    static bool valid_input(boostPO::variables_map&);

private:

    void get_tax_scope(std::string&, QuerySequence::EggnogResults&);
    void get_sql_data(QuerySequence::EggnogResults&, SQLDatabaseHelper&);
    std::string format_sql_data(std::string&);
    void get_og_query(QuerySequence::EggnogResults&);

    static constexpr short EGGNOG_COL_NUM     = 12;
    static constexpr uint8 COUNT_TOP_TAX_SCOPE= 10;
    static constexpr uint8 COUNT_TOP_GO       = 10;
    const std::string EGGNOG_DIRECTORY        = "EggNOG/";
    const std::string GRAPH_EGG_TAX_BAR_TITLE = "Top_10_Tax_Levels";
    const std::string GRAPH_EGG_TAX_BAR_PNG   = "eggnog_tax_scope.png";
    const std::string GRAPH_EGG_TAX_BAR_TXT   = "eggnog_tax_scope.txt";
    const std::string EGG_ANNOT_RESULTS       = "annotation_results";
    const std::string EGG_ANNOT_STD           = "annotation_std";
    const std::string EGG_ANNOT_APPEND        = ".emapper.annotations";

    std::string mFigureDir;
    std::string mProcDir;
    std::string _egg_out_dir;
    std::string mEggnogDbPath;
    std::string mOutHIts;

    static std::string EGG_EMAPPER_EXE;
    std::string eggnog_format(std::string);
};


#endif //ENTAP_MODEGGNOG_H
#endif