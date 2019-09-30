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


//*********************** Includes *****************************
#include "EntapGlobals.h"
//**************************************************************

// This table should match EntapGlobals.h ENTAP_HEADERS enum
EntapHeader ENTAP_HEADER_INFO[] = {
        {"Unused",              false},                         // 0
        {"Query Sequence",      true},

        /* Frame Selection */
        {"Frame",               true},

        /* Expression Filtering */
        {"FPKM",                true},
        {"TPM",                 true},

        /* Similarity Search - General */
        {"Subject Sequence",    true},
        {"Percent Identical",   true},                          // 5
        {"Alignment Length",    true},
        {"Mismatches",          true},
        {"Gap Openings",        true},
        {"Query Start",         true},
        {"Query End",           true},                          // 10
        {"Subject Start",       true},
        {"Subject End",         true},
        {"E Value",             true},
        {"Coverage",            true},
        {"Description",         true},                          // 15
        {"Species",             true},
        {"Taxonomic Lineage",   true},
        {"Origin Database",     true},
        {"Contaminant",         true},
        {"Informative",         true},

        /* Similarity Search - UniProt */
        {"UniProt Database Cross Reference",        true},      // 20
        {"UniProt Additional Information",          true},
        {"UniProt KEGG Terms",                      true},
        {"UniProt GO Biological",                   true},
        {"UniProt GO Cellular",                     true},
        {"UniProt GO Molecular",                    true},

        /* Ontology - EggNOG */
        {"EggNOG Seed Ortholog",                    true},
        {"EggNOG Seed E-Value",                     true},
        {"EggNOG Seed Score",                       true},
        {"EggNOG Predicted Gene",                   true},
        {"EggNOG Tax Scope",                        true},
        {"EggNOG Tax Scope Max",                    true},
        {"EggNOG Member OGs",                       true},
        {"EggNOG EggNOG Description",               false},
        {"EggNOG BIGG Reaction",                    true},
        {"EggNOG KEGG Terms",                       true},
        {"EggNOG GO Biological",                    true},
        {"EggNOG GO Cellular",                      true},
        {"EggNOG GO Molecular" ,                    true},
        {"EggNOG Protein Domains",                  false},

        /* Ontology - InterProScan */
        {"IPScan GO Biological",                    true},
        {"IPScan GO Cellular",                      true},
        {"IPScan GO Molecular",                     true},
        {"IPScan Pathways",                         true},
        {"IPScan InterPro ID",                      true},
        {"IPScan Protein Database",                 true},
        {"IPScan Protein Description",              true},
        {"IPScan E-Value",                          true},


        {"Unused",                                  false}
};

// For C++11 use std::to_string(arg)
std::string float_to_string(fp64 val) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << val;
    return ss.str();
}

std::string float_to_sci(fp64 val, int precision) {
    std::stringstream ss;
    ss << std::fixed <<
          std::setprecision(precision) <<
          std::scientific              <<
          val;
    return ss.str();
}

vect_str_t split_string(std::string sequences, char delim) {
    vect_str_t split_vals;

    // Remove newline if exists
    sequences.erase(std::remove(sequences.begin(), sequences.end(), '\n'), sequences.end());

    std::istringstream iss(sequences);
    std::string val;
    while(std::getline(iss, val, delim)) {
        split_vals.push_back(val);
    }
    return split_vals;
}

std::string get_cur_time() {
    std::chrono::time_point<std::chrono::system_clock> current;
    std::time_t time;

    current = std::chrono::system_clock::now();
    time = std::chrono::system_clock::to_time_t(current);
    std::string out_time(std::ctime(&time));
    return out_time.substr(0,out_time.length()-1);
}