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

#ifndef ENTAP_ENTAPMODULE_H
#define ENTAP_ENTAPMODULE_H

#include "common.h"
#include "EntapGlobals.h"
#include "UserInput.h"
#include "database/EntapDatabase.h"

class EntapModule {

public:

    struct ModVerifyData {
        bool        files_exist;
        vect_str_t  output_paths;
    };


    EntapModule(std::string &execution_stage_path, std::string &in_hits,
                EntapDataPtrs &entap_data, std::string module_name, std::string &exe);

    virtual ~EntapModule() = default;
    virtual ModVerifyData verify_files()=0;
    virtual void execute() = 0;
    virtual void parse() = 0;

protected:

    const std::string PROCESSED_OUT_DIR     = "processed/";
    const std::string FIGURE_DIR            = "figures/";
    const std::string OVERALL_RESULTS_DIR   = "overall_results";

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

    const std::string YES_FLAG                = "Yes";
    const std::string NO_FLAG                 = "No";

    const uint8 COUNT_TOP_GO                  = 10;
    const uint8 COUNT_TOP_SPECIES             = 10;

    bool               _blastp;
    bool               _overwrite;
    int                _threads;
    uint16             _software_flag;
    std::string        _outpath;
    std::string        _in_hits;
    std::string        _proc_dir;                   // "processed" directory, or data analyzed
    std::string        _figure_dir;
    std::string        _mod_out_dir;
    std::string        _overall_results_dir;
    std::string        _exe_path;
    std::string        _transcript_shortname;       // filename of transcriptome file for file name purposes
    ExecuteStates      _execution_state;
    std::vector<uint16> _go_levels;
    GraphingManager    *_pGraphingManager;
    QueryData          *_pQUERY_DATA;
    UserInput          *_pUserInput;
    FileSystem         *_pFileSystem;
    EntapDatabase      *_pEntapDatabase;

    go_format_t EM_parse_go_list(std::string list, EntapDatabase* database,char delim);
    void EM_set_uniprot_headers(bool set);
};


#endif //ENTAP_ENTAPMODULE_H
