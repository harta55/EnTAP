/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
<<<<<<< HEAD
 * Copyright 2017-2020, Alexander Hart, Dr. Jill Wegrzyn
=======
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
>>>>>>> master
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
<<<<<<< HEAD
                EntapDataPtrs &entap_data, std::string module_name,
                std::vector<ENTAP_HEADERS> &module_headers);
=======
                EntapDataPtrs &entap_data, std::string module_name, std::string &exe);

>>>>>>> master
    virtual ~EntapModule() = default;
    virtual ModVerifyData verify_files()=0;
    virtual void execute() = 0;
    virtual void parse() = 0;
<<<<<<< HEAD
    virtual void set_success_flags() = 0;
    virtual void get_version() = 0;
=======
>>>>>>> master

protected:

    const std::string PROCESSED_OUT_DIR     = "processed/";
    const std::string FIGURE_DIR            = "figures/";
    const std::string OVERALL_RESULTS_DIR   = "overall_results";

    const std::string FILENAME_OUT_UNANNOTATED  = "unannotated_sequences";
    const std::string FILENAME_OUT_ANNOTATED    = "annotated_sequences";

    const std::string GRAPH_GO_END_TXT        = "_go_bar_graph.txt";
    const std::string GRAPH_GO_END_PNG        = "_go_bar_graph.png";
    const std::string GRAPH_GO_BAR_BIO_TITLE  = "Top_GO_Biological_Terms";
    const std::string GRAPH_GO_BAR_CELL_TITLE = "Top_GO_Cellular_Terms";
    const std::string GRAPH_GO_BAR_MOLE_TITLE = "Top_GO_Molecular_Terms";
    const std::string GRAPH_GO_BAR_ALL_TITLE  = "Top_GO_Terms";

    const std::string YES_FLAG                = "Yes";
    const std::string NO_FLAG                 = "No";

    const uint16 COUNT_TOP_GO                  = 10;
    const uint16 COUNT_TOP_SPECIES             = 10;

<<<<<<< HEAD
    static std::vector<ENTAP_HEADERS>   mModuleHeaders;
    bool               mBlastp;                 // TRUE/FALSE whether user has input blastp
    bool               mOverwrite;              // Indicates whether user would like to delete previous execution files for module
    int                mThreads;                // Number of threads specified for execution from user
    uint16             mSoftwareFlag;           // Flag indicating software module being used
    uint16             mVersionMajor;
    uint16             mVersionMinor;
    uint16             mVersionRev;
    std::string        mOutpath;
    ent_input_str_t    mInputTranscriptome;
    std::string        mModuleName;                // Name of module
    std::string        mProcDir;                   // "processed" directory, or data analyzed
    std::string        mFigureDir;                 // "figure" directory to place any generated figures
    std::string        mModOutDir;                 // Root out directory for the EnTAP module we are dealing with
    std::string        mOverallResultsDir;
    ent_input_str_t    mExePath;
    std::string        mTranscriptomeShortname;       // filename of transcriptome file for file name purposes
    ExecuteStates      mExecutionState;
    std::vector<uint16> mGoLevels;
    std::vector<ENTAP_HEADERS> mEntapHeaders;
    GraphingManager    *mpGraphingManager;
    QueryData          *mpQueryData;
    UserInput          *mpUserInput;
    FileSystem         *mpFileSystem;
    EntapDatabase      *mpEntapDatabase;
    std::vector<FileSystem::ENT_FILE_TYPES> mAlignmentFileTypes; // may be overriden by module

    go_format_t EM_parse_go_list(std::string list, EntapDatabase* database,char delim);
    void enable_headers();
    virtual std::vector<ENTAP_HEADERS>  &moduleHeaders() { return mModuleHeaders; }
=======
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
    std::vector<FileSystem::ENT_FILE_TYPES> _alignment_file_types; // may be overriden by module

    go_format_t EM_parse_go_list(std::string list, EntapDatabase* database,char delim);
>>>>>>> master
};


#endif //ENTAP_ENTAPMODULE_H
