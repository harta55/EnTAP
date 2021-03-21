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

#ifndef ENTAP_MODBUSCO_H
#define ENTAP_MODBUSCO_H


#include "AbstractOntology.h"

/*
 * BUSCO Reference:
 *
 * Seppey M., Manni M., Zdobnov E.M. (2019) BUSCO: Assessing Genome Assembly and
 * Annotation Completeness. In: Kollmar M. (eds) Gene Prediction. Methods in Molecular
 * Biology, vol 1962. Humana, New York, NY. 2019 doi.org/10.1007/978-1-4939-9173-0_14.
 * PMID:31020564
 */

class ModBUSCO : public AbstractOntology {
public:
    ModBUSCO(std::string &in_hits, std::string &ont_out, EntapDataPtrs &entap_data);
    ~ModBUSCO();
    virtual void execute() override;
    virtual void parse() override;
    virtual ModVerifyData verify_files() override;
    static bool is_executable(std::string &exe);
    virtual bool set_version() override ;

private:


    typedef enum {
        BUSCO_VERSION_UNKNOWN,      // Not supported
        BUSCO_VERSION_3,            // Not supported
        BUSCO_VERSION_4             // Supported :)
    } BUSCO_VERSION;

    std::string get_final_table_path(std::string &output_dir, std::string &database, BUSCO_VERSION &version);
    void set_version_defaults();
    bool is_version_supported(BUSCO_VERSION &version);
    std::string print_version();

    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;

    // WARNING commands current as of version 3.0.2 and 4.0.2, subject to change!!
    const std::string BUSCO_INPUT_IN       = "--in";        // Specify input transcriptome
    const std::string BUSCO_INPUT_OUTPUT   = "-o";          // Specify output flag/directory (sort of works?)
    const std::string BUSCO_INPUT_RUN_TYPE = "-m";          // Specify execution type (prot, tran)
    const std::string BUSCO_INPUT_DATABASE = "-l";          // Specify path to BUSCO database
    const std::string BUSCO_INPUT_CPU      = "--cpu";       // Specify number of threads to BUSCO
    const std::string BUSCO_INPUT_EVAL     = "--evalue";

    // BUSCO 4.0.2 commands
    const std::string BUSCO_RUN_TYPE_TRAN  = "transcriptome";
    const std::string BUSCO_RUN_TYPE_PROT  = "proteins";

    // BUSCO 3.0.2 prepends this tag before all directory creations (-o command)
    // BUSCO 4.0.2 does not do this
    const std::string BUSCO_PREPEND_TAG    = "run_";
    static constexpr uint16 BUSCO_COLUMN_NUM     = 5;         // full_table.tsv column number
    const std::string BUSCO_FULL_TABLE_FILENAME = "full_table.tsv";
    const std::string BUSCO_STATUS_COMPLETE= "Complete";
    const std::string BUSCO_STATUS_MISSING = "Missing";

    const BUSCO_VERSION              DEFAULT_VERSION    = BUSCO_VERSION_4;
    const std::vector<BUSCO_VERSION> SUPPORTED_VERSIONS = {BUSCO_VERSION_4};
    const uint16                     DEFAULT_VERSION_MAJOR = 4;
    const uint16                     DEFAULT_VERSION_MINOR = 0;
    const uint16                     DEFAULT_VERSION_REV = 2;

    BUSCO_VERSION    mBuscoVersion;
    ent_input_str_t  mExePath;       // Path to BUSCO executable
    std::string      mRunType;       // Execution type for BUSCO (prot or tran based on blastp)
    std::string      mFinalTablePath;
    std::string      mOutputDirectoryTag; // output directory TAG that will go into BUSCO as -o command for particular database
    std::string      mOutputRunDir;    // Full path to output directory for busco run with TAG
    ent_input_str_t  mBuscoDatabase; // BUSCO database we want to use
    ent_input_fp_t   mEval;          // Eval to use during BUSCO search


};


#endif //ENTAP_MODBUSCO_H
