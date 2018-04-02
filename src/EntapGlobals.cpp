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


//*********************** Includes *****************************
#include <chrono>
#include <ios>
#include <boost/filesystem/operations.hpp>
#include <pstream.h>
#include <boost/program_options/variables_map.hpp>
#include <unordered_map>
#include <boost/archive/binary_iarchive.hpp>
#include "EntapGlobals.h"
#include "config.h"
#include "FileSystem.h"
#include <boost/archive/binary_iarchive.hpp>
#include <iomanip>

//**************************************************************


namespace ENTAP_EXECUTE {
    //------------------------Ontology-------------------------//
    const std::string GO_BIOLOGICAL_FLAG = "biological_process";
    const std::string GO_CELLULAR_FLAG   = "cellular_component";
    const std::string GO_MOLECULAR_FLAG  = "molecular_function";


    //------------------------Headers-------------------------//

    // Similarity Search
    const std::string HEADER_QUERY     = "Query Seq";
    const std::string HEADER_SUBJECT   = "Subject Seq";
    const std::string HEADER_PERCENT   = "Percent Identical";
    const std::string HEADER_ALIGN_LEN = "Alignment Length";
    const std::string HEADER_MISMATCH  = "Mismatches";
    const std::string HEADER_GAP_OPEN  = "Gap Openings";
    const std::string HEADER_QUERY_S   = "Query Start";
    const std::string HEADER_QUERY_E   = "Query End";
    const std::string HEADER_SUBJ_S    = "Subject Start";
    const std::string HEADER_SUBJ_E    = "Subject End";
    const std::string HEADER_E_VAL     = "E Value";
    const std::string HEADER_COVERAGE  = "Coverage";
    const std::string HEADER_TITLE     = "Description";
    const std::string HEADER_SPECIES   = "Species";
    const std::string HEADER_DATABASE  = "Origin Database";
    const std::string HEADER_FRAME     = "Frame";
    const std::string HEADER_CONTAM    = "Contaminant";
    const std::string HEADER_INFORM    = "Informative";

    // EggNOG
    const std::string HEADER_SEED_ORTH   = "Seed Ortholog";
    const std::string HEADER_SEED_EVAL   = "Seed E-Value";
    const std::string HEADER_SEED_SCORE  = "Seed Score";
    const std::string HEADER_PRED_GENE   = "Predicted Gene";
    const std::string HEADER_TAX_SCOPE   = "Tax Scope";
    const std::string HEADER_EGG_OGS     = "OGs";
    const std::string HEADER_EGG_KEGG    = "KEGG Terms";
    const std::string HEADER_EGG_GO_BIO  = "GO Biological";
    const std::string HEADER_EGG_GO_CELL = "GO Cellular";
    const std::string HEADER_EGG_GO_MOLE = "GO Molecular";
    const std::string HEADER_EGG_DESC    = "Eggnog Description";
    const std::string HEADER_EGG_LEVEL   = "Full Tax Scope";
    const std::string HEADER_EGG_PROTEIN = "Protein Domains";

    // InterProScan
    const std::string HEADER_INTER_GO_BIO    = "IPScan GO Biological";
    const std::string HEADER_INTER_GO_CELL   = "IPScan GO Cellular";
    const std::string HEADER_INTER_GO_MOLE   = "IPScan GO Molecular";
    const std::string HEADER_INTER_PATHWAY   = "Pathways";
    const std::string HEADER_INTER_INTERPRO  = "InterPro";
    const std::string HEADER_INTER_DATA_TYPE = "Protein Database";
    const std::string HEADER_INTER_DATA_TERM = "Protein Description";
    const std::string HEADER_INTER_EVAL      = "E Value";
}

namespace UInput {
    //------------------USER INPUTS-----------------------//
    const std::string INPUT_FLAG_TAG           = "out-dir";
    const std::string INPUT_FLAG_CONFIG        = "config";
    const std::string INPUT_FLAG_ALIGN         = "align";
    const std::string INPUT_FLAG_RUNPROTEIN    = "runP";
    const std::string INPUT_FLAG_RUNNUCLEOTIDE = "runN";
    const std::string INPUT_FLAG_OVERWRITE     = "overwrite";
    const std::string INPUT_FLAG_NCBI_1        = "ncbi";
    const std::string INPUT_FLAG_NCBI_2        = "N";
    const std::string INPUT_FLAG_UNIPROT       = "uniprot";
    const std::string INPUT_FLAG_INTERPRO      = "protein";
    const std::string INPUT_FLAG_ONTOLOGY      = "ontology";
    const std::string INPUT_FLAG_SPECIES       = "taxon";
    const std::string INPUT_FLAG_QCOVERAGE     = "qcoverage";
    const std::string INPUT_FLAG_TCOVERAGE     = "tcoverage";
    const std::string INPUT_FLAG_COMPLETE      = "complete";
    const std::string INPUT_FLAG_GO_LEVELS     = "level";
    const std::string INPUT_FLAG_EXE_PATH      = "paths";
    const std::string INPUT_FLAG_FPKM          = "fpkm";
    const std::string INPUT_FLAG_DATA_OUT      = "database-out";
    const std::string INPUT_FLAG_CONTAM        = "contam";
    const std::string INPUT_FLAG_E_VAL         = "e";
    const std::string INPUT_FLAG_HELP          = "help";
    const std::string INPUT_FLAG_VERSION       = "version";
    const std::string INPUT_FLAG_TRANSCRIPTOME = "input";
    const std::string INPUT_FLAG_DATABASE      = "database";
    const std::string INPUT_FLAG_GRAPH         = "graph";
    const std::string INPUT_FLAG_TRIM          = "trim";
    const std::string INPUT_FLAG_STATE         = "state";
    const std::string INPUT_FLAG_SINGLE_END    = "single-end";
    const std::string INPUT_FLAG_THREADS       = "threads";
    const std::string INPUT_FLAG_UNINFORM      = "uninformative";
    const std::string INPUT_FLAG_NOCHECK       = "no-check";

}

/**
 * ======================================================================
 * Function int execute_cmd(std::string cmd, std::string out_path)
 *
 * Description          - Terminal stream based on pstreams implementation
 *                      - Executes commands and prints .err and .out from stream
 *
 * Notes                - None
 *
 * @param cmd           - Command for child process
 * @param out_path      - Path to std out/err files to be printed
 *
 * @return              - int error code
 *
 * =====================================================================
 */
int execute_cmd(std::string cmd, std::string out_path) {
    std::ofstream out_file(out_path+".out", std::ios::out | std::ios::app);
    std::ofstream err_file(out_path+".err", std::ios::out | std::ios::app);
    const redi::pstreams::pmode mode = redi::pstreams::pstdout|redi::pstreams::pstderr;
    redi::ipstream child(cmd, mode);
    char buf[1024];
    std::streamsize n;
    bool finished[2] = { false, false };
    while (!finished[0] || !finished[1]) {
        if (!finished[0]) {
            while ((n = child.err().readsome(buf, sizeof(buf))) > 0)
                err_file.write(buf, n);
            if (child.eof()) {
                finished[0] = true;
                if (!finished[1])
                    child.clear();
            }
        }
        if (!finished[1]) {
            while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                out_file.write(buf, n).flush();
            if (child.eof()) {
                finished[1] = true;
                if (!finished[0])
                    child.clear();
            }
        }
    }
    child.close();
    out_file.close();
    err_file.close();
    if (child.rdbuf()->exited())
        return child.rdbuf()->status();
    return 1;
}
// todo, may want to handle differently
// TODO change to sending map of flags as command
int execute_cmd(std::string cmd) {
    const redi::pstreams::pmode mode = redi::pstreams::pstdout|redi::pstreams::pstderr;
    redi::ipstream child(cmd, mode);
    char buf[1024];
    std::streamsize n;
    bool finished[2] = { false, false };
    while (!finished[0] || !finished[1]) {
        if (!finished[0]) {
            while ((n = child.err().readsome(buf, sizeof(buf))) > 0)
                continue;
            if (child.eof()) {
                finished[0] = true;
                if (!finished[1])
                    child.clear();
            }
        }
        if (!finished[1]) {
            while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                continue;
            if (child.eof()) {
                finished[1] = true;
                if (!finished[0])
                    child.clear();
            }
        }
    }
    child.close();
    if (child.rdbuf()->exited())
        return child.rdbuf()->status();
    return 1;
}


std::string generate_command(std::unordered_map<std::string,std::string> &map,std::string exe_path) {
    std::stringstream ss;
    std::string       out;

    ss << exe_path << " ";
    for (auto &pair : map)ss << pair.first << " " << pair.second << " ";
    out = ss.str();
    return out;
}

std::string float_to_string(fp64 val) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << val;
    return ss.str();
}
