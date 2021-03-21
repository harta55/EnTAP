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


#ifndef ENTAP_MODINTERPRO_H
#define ENTAP_MODINTERPRO_H


#include "AbstractOntology.h"
#include "../UserInput.h"
#include "../QuerySequence.h"

/**
 * ======================================================================
 * @class ModInterpro
 *
 * Description          - This EnTAP module supports execution, parsing, and
 *                        statistical analysis of the InterProScan software
 *                        through terminal commands
 *                      - Parsed data is added to QueryData class
 *                      - Inherits from AbstractOntology and EntapModule classes
 *
 * Citation             - P. Jones et al., “InterProScan 5: genome-scale
 *                        protein function classification,” (in eng),
 *                        Bioinformatics, vol. 30, no. 9, pp. 1236-40, May 2014.
 *
 * ======================================================================
 */
class ModInterpro : public AbstractOntology{

    struct InterProData {
        std::string interID;
        std::string interDesc;
        std::string databaseID;
        std::string databasetype;
        std::string databaseDesc;
        std::string pathways;
        std::string go_terms;
        fp64        eval;
    };

public:
    ~ModInterpro();
    ModInterpro(std::string &ont, std::string &in, EntapDataPtrs& entap_data);

    virtual ModVerifyData verify_files() override ;
    virtual void execute() override ;
    virtual void parse() override ;
    static bool is_executable(std::string &exe);
    virtual bool set_version() override;

    static bool valid_input(UserInput* userinput);
private:
    std::string NUCLEO_TAG                  = "n";
    std::string PROTEIN_TAG                 = "p";
    std::string INTERPRO_TEMP               = "temp/";
    std::string INTERPRO_OUTPUT             = "interpro_results";
    std::string INTERPRO_STD_OUT            = "interproscan_std";
    std::string OUT_HITS_FAA                = "interpro_hits.faa";
    std::string OUT_HITS_FNN                = "interpro_hits.fnn";
    std::string OUT_NO_HITS_FAA             = "interpro_no_hits.faa";
    std::string OUT_NO_HITS_FNN             = "interpro_no_hits.fnn";
    std::string INTERPRO_EXT_XML            = ".xml";
    std::string INTERPRO_EXT_TSV            = ".tsv";

    std::string INTERPRO_DATABASE_FLAG              = "interpro";   // TODO add full alignment support
    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;

    // Valid databases
    static const std::vector<std::string> INTERPRO_DATABASES;
    static const std::string INTERPRO_DEFAULT;

    static constexpr short INTERPRO_COL_NUM = 15;

    std::string XML_SIGNATURE = "signature";
    std::string XML_ENTRY     = "entry";
    std::string XML_XREF      = "xref";
    std::string XML_PRO_M     = "protein-matches";
    std::string XML_PROTEIN   = "protein";
    std::string XML_MATCHES   = "matches";

    std::string FLAG_SEQTYPE  = " --seqtype";
    std::string FLAG_GOTERM   = " --goterms";
    std::string FLAG_IPRLOOK  = " --iprlookup";
    std::string FLAG_PATHWAY  = " --pathways";
    std::string FLAG_TEMP     = " --tempdir";

    vect_str_t               mDatabases;
    std::string              mFinalOutpath;
    std::string              mFinalBasepath;

#if 0
    std::map<std::string,InterProData> parse_xml(void);
#endif
    std::map<std::string,InterProData> parse_tsv();
    std::string format_interpro();
};


#endif //ENTAP_MODINTERPRO_H
