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


#ifndef ENTAP_ABSTRACTONTOLOGY_H
#define ENTAP_ABSTRACTONTOLOGY_H

#include <string>
#include "../GraphingManager.h"
#include "../QuerySequence.h"
#include "../QueryData.h"
#include "../database/EntapDatabase.h"
#include "../EntapGlobals.h"

class QuerySequence;
class QueryData;
class GoTerm;
class AbstractOntology {
public:
    AbstractOntology(std::string &out,
                     std::string &in_hits,
                     std::string &ont_out,
                     bool blastp,
                     std::vector<uint16>& lvls,
                     EntapDataPtrs &entap_data){
        _outpath         = out;
        _inpath          = in_hits;
        _ontology_dir    = ont_out;
        _blastp          = blastp;
        _go_levels       = lvls;

        _pGraphingManager = entap_data._pGraphingManager;
        _pQUERY_DATA      = entap_data._pQueryData;
        _pFileSystem      = entap_data._pFileSystem;
        _pUserInput       = entap_data._pUserInput;
        _pEntapDatabase   = entap_data._pEntapDatbase;

        _threads         = _pUserInput->get_supported_threads();
    }

    virtual ~AbstractOntology() = default;
    virtual std::pair<bool, std::string> verify_files()=0;
    virtual void execute() = 0;
    virtual void parse() = 0;

protected:
    const std::string PROCESSED_OUT_DIR     = "processed/";
    const std::string FIGURE_DIR            = "figures/";

    const std::string GO_MOLECULAR_FLAG     = "molecular_function";
    const std::string GO_BIOLOGICAL_FLAG    = "biological_process";
    const std::string GO_CELLULAR_FLAG      = "cellular_component";
    const std::string GO_OVERALL_FLAG       = "overall";

    const std::string OUT_UNANNOTATED_NUCL  = "unannotated_sequences.fnn";
    const std::string OUT_UNANNOTATED_PROT  = "unannotated_sequences.faa";
    const std::string OUT_ANNOTATED_NUCL    = "annotated_sequences.fnn";
    const std::string OUT_ANNOTATED_PROT    = "annotated_sequences.faa";

    const std::string GRAPH_GO_END_TXT        = "_go_bar_graph.txt";
    const std::string GRAPH_GO_END_PNG        = "_go_bar_graph.png";
    const std::string GRAPH_GO_BAR_BIO_TITLE  = "Top_10_GO_Biological_Terms";
    const std::string GRAPH_GO_BAR_CELL_TITLE = "Top_10_GO_Cellular_Terms";
    const std::string GRAPH_GO_BAR_MOLE_TITLE = "Top_10_GO_Molecular_Terms";
    const std::string GRAPH_GO_BAR_ALL_TITLE  = "Top_10_GO_Terms";

    const uint8 GRAPH_ONTOLOGY_FLAG = 4;
    const uint8 GRAPH_TOP_BAR_FLAG  = 1;  // used for tax levels and go term tops

    bool               _blastp;
    int                _threads;
    std::string        _outpath;
    std::string        _inpath;
    std::string        _ontology_dir;
    std::vector<uint16> _go_levels;
    GraphingManager    *_pGraphingManager;
    QueryData          *_pQUERY_DATA;
    UserInput          *_pUserInput;
    FileSystem         *_pFileSystem;
    EntapDatabase      *_pEntapDatabase;

    std::map<std::string,std::vector<std::string>> parse_go_list(std::string list, EntapDatabase*,char delim);
};


#endif //ENTAP_ABSTRACTONTOLOGY_H
