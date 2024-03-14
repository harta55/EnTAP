/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2024, Alexander Hart, Dr. Jill Wegrzyn
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
#include "../QuerySequence.h"
#include "../QueryData.h"
#include "../GraphingManager.h"
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
    ModEggnog(std::string &ont_out, std::string &in_hits,
            EntapDataPtrs &entap_data);
    ModVerifyData verify_files() override;
    ~ModEggnog();
    virtual void execute() override ;
    virtual void parse() override;
    static bool is_executable(std::string &exe);
    virtual bool set_version() override;

private:

    enum EGGNOG_MAPPER_STATES {
        EGGNOG_MAPPER_NOT_STARTED=0,
        EGGNOG_MAPPER_DIAMOND_COMPLETE,
        EGGNOG_MAPPER_ANNOTATIONS_COMPLETE
    };


    static constexpr short EGGNOG_COL_NUM     = 21;
    static constexpr uint8 COUNT_TOP_TAX_SCOPE= 10;
    static constexpr uint8 COUNT_TOP_GO       = 10;
    const std::string EGGNOG_DIRECTORY        = "EggNOG/";
    const std::string GRAPH_EGG_TAX_BAR_TITLE = "Top_10_Tax_Levels";
    const std::string GRAPH_EGG_TAX_BAR_PNG   = "eggnog_tax_scope.png";
    const std::string GRAPH_EGG_TAX_BAR_TXT   = "eggnog_tax_scope.txt";
    const std::string EGG_ANNOT_RESULTS       = "annotation_results";
    const std::string EGG_ANNOT_STD           = "annotation_std";
    const std::string EGG_OUTPUT_ANNOT_APPEND = ".emapper.annotations";
    const std::string EGG_OUTPUT_HITS_APPEND  = ".emapper.hits";
    const std::string EGG_OUTPUT_SEED_ORTHO_APPEND = ".emapper.seed_orthologs";
    const std::string EGG_OUT_UNANNOTATED     = "unannotated";
    const std::string EGG_OUT_ANNOTATED       = "annotated";
    const std::string EGG_MAPPER_PREFIX   = "eggnog_";

    std::string mOutHIts;
    std::string mEggnogMapDataDir;
    std::string mEggnogMapDMNDPath;
    std::string mEggnogMapAnnotationsOutputPath;
    std::string mEggnogMapHitsOutputPath;
    std::string mEggnogMapSeedOrthoOutputPath;
    EGGNOG_MAPPER_STATES mEggnogMapperState;

    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;
    std::string get_output_tag();
};


#endif //ENTAP_MODEGGNOG_H