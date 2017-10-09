/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
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
#include "ModInterpro.h"
#include "../ExceptionHandler.h"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/xml_parser.hpp"
//**************************************************************

using boost::property_tree::ptree;


/**
 * ======================================================================
 * Function std::pair<bool, std::string> ModInterpro::verify_files()
 *
 * Description          - Checks whether execution has already been ran
 *                        with the same input (so it can be skipped)
 *
 * Notes                - None
 *
 *
 * @return              - Pair of yes/no if files were found and string (not used)
 *
 * =====================================================================
 */
std::pair<bool, std::string> ModInterpro::verify_files() {

    std::string filename;

    filename       = INTERPRO_OUTPUT;
    _final_basepath= PATHS(_interpro_dir, filename);
    filename      += INTERPRO_EXT;
    _final_outpath = PATHS(_interpro_dir, filename);
    return std::make_pair(file_exists(_final_outpath), "");
}

void ModInterpro::execute() {

    std::string interpro_cmd;
    std::string std_out;

    std_out      = PATHS(_interpro_dir, INTERPRO_STD_OUT);
    interpro_cmd =
            INTERPRO_EXE    +
            " -i "          + _inpath +
            " -b "          + _final_basepath +
            FLAG_GOTERM     +
            FLAG_IPRLOOK    +
            FLAG_PATHWAY;

    print_debug("Executing InterProScan:\n" + interpro_cmd);
    if (!_databases.empty()) {
        for (std::string &val : _databases) interpro_cmd += " --appl " + val;
    } else {
        throw ExceptionHandler("No InterPro databases selected!",
                ENTAP_ERR::E_RUN_INTERPRO);
    }
    if (execute_cmd(interpro_cmd, std_out) != 0) {
        throw ExceptionHandler("Error executing InterProScan, consult the error file at: "+
                std_out, ENTAP_ERR::E_RUN_INTERPRO);
    }
}


/**
 * ======================================================================
 * Function void ModInterpro::parse()
 *
 * Description          - Parses XML file produced from InterProScan
 *                      - Will probably be changed with new library introduction
 *
 * Notes                - None
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModInterpro::parse() {

    std::stringstream                     stats_stream;
    std::string                           stats_out;
    std::string                           e_str;
    std::string                           interpro_output;
    std::string                           protein_output;
    std::map<std::string, struct_go_term> GO_DATABASE;
    std::map<std::string,InterProData>    interpro_map;
    go_struct                             go_terms_parsed;
    std::string                           seq_id;
    fp64                                  e_val;
    bool                                  inter;
    std::string                           pathway;
    ptree                                 pt;

    print_debug("Beginning to parse InterProScan data...");

    try {
        GO_DATABASE = read_go_map();
    } catch (ExceptionHandler const &e) {throw e;}

    if (file_exists(_final_outpath)) {
        boost::property_tree::read_xml(_final_outpath, pt);
        print_debug("File found at: " + _final_outpath + " parsing...");
    } else {
        throw ExceptionHandler("Unable to locate InterProScan file at: " +
            _final_outpath, ENTAP_ERR::E_PARSE_INTERPRO);
    }
    for (ptree::value_type const& v : pt.get_child(XML_PRO_M)) {
        if (v.first == XML_PROTEIN) {
            seq_id = v.second.get_child(XML_XREF).get("<xmlattr>.id","");
            for (ptree::value_type const& s : v.second.get_child(XML_MATCHES)){
                if (s.first.find("-match") != std::string::npos) {
                    e_val = s.second.get("<xmlattr>.evalue", 1.0);
                    if (interpro_map.find(seq_id) != interpro_map.end()) {
                        if (interpro_map[seq_id].eval < e_val) continue;
                    }
                    InterProData interpro1;
                    interpro1.pathways = "";
                    interpro1.eval     = e_val;
                    try {
                        s.second.get_child(XML_SIGNATURE).get_child(XML_ENTRY);
                        inter = true;
                    } catch (...) {inter = false;}
                    interpro1.databaseDesc = s.second.get_child(XML_SIGNATURE).get("<xmlattr>.desc","");
                    interpro1.databaseID   = s.second.get_child(XML_SIGNATURE).get("<xmlattr>.ac","");
                    interpro1.databasetype = s.second.get_child(XML_SIGNATURE).get_child("signature-library-release").
                            get("<xmlattr>.library", "");
                    if (inter) {
                        interpro1.interDesc = s.second.get_child(XML_SIGNATURE).get_child(XML_ENTRY).
                                get("<xmlattr>.desc", "");
                        interpro1.interID = s.second.get_child(XML_SIGNATURE).get_child(XML_ENTRY).
                                get("<xmlattr>.ac", "");
                        std::unordered_map<std::string,std::string> pathway_map;
                        for (ptree::value_type const& t : s.second.get_child(XML_SIGNATURE).get_child(XML_ENTRY)) {
                            if (t.first.find("go-xref") == 0) {
                                interpro1.go_terms += t.second.get("<xmlattr>.id","") + ",";
                            }
                            if (interpro1.go_terms.length() > 1) interpro1.go_terms.pop_back();
                            if (t.first.find("pathway-xref") == 0) {
                                pathway = t.second.get("<xmlattr>.db","");
                                if (pathway_map.find(pathway)!=pathway_map.end()) {
                                    pathway_map[pathway] += ", ";
                                }
                                pathway_map[pathway] +=
                                        t.second.get("<xmlattr>.id","")      +  "-"  +
                                        t.second.get("<xmlattr>.name","");
                            }
                        }
                        for (auto &entry : pathway_map) {
                            if (entry.first.empty()) continue;
                            interpro1.pathways += entry.first + "(" + entry.second + ");";
                        }
                        if (interpro1.pathways.length() > 1) interpro1.pathways.pop_back();
                    } else {
                        interpro1.interDesc = "";
                        interpro1.interID = "";
                    }
                    interpro_map[seq_id] = interpro1;
                }
            }
        }
    }

    print_debug("Success! Beginning to update query sequences...");

    // TODO stats
    for (auto &pair : *pQUERY_DATA->get_sequences_ptr()) {
        std::map<std::string, InterProData>::iterator it = interpro_map.find(pair.first);
        if (it != interpro_map.end()) {
            interpro_output = it->second.interID + "(" + it->second.interDesc + ")";
            protein_output  = it->second.databaseID + "(" + it->second.databaseDesc + ")";
            go_terms_parsed = parse_go_list(it->second.go_terms,GO_DATABASE,',');
            std::stringstream ss;
            ss << std::scientific << it->second.eval;
            e_str = ss.str();
            pair.second.set_interpro_results(e_str, protein_output, it->second.databasetype,
                                             interpro_output, it->second.pathways, go_terms_parsed);
        } else {
            ;
        }
    }
    print_debug("Success! Calculating statistics...");
    stats_stream <<
                 ENTAP_STATS::SOFTWARE_BREAK << " Ontology - InterProScan" << ENTAP_STATS::SOFTWARE_BREAK <<
                 "InterProScan statistics coming soon!";
    stats_out = stats_stream.str();
    print_statistics(stats_out);
    print_debug("Success! InterProScan finished");
}

void ModInterpro::set_data(std::string & unused, std::vector<std::string>& interpro) {
    _interpro_dir = PATHS(_ontology_dir, INTERPRO_DIRECTORY);
    _proc_dir     = PATHS(_interpro_dir, PROCESSED_OUT_DIR);
    _figure_dir   = PATHS(_interpro_dir, FIGURE_DIR);
    _databases    = interpro;

    boostFS::remove_all(_proc_dir);
    boostFS::remove_all(_figure_dir);

    boostFS::create_directories(_interpro_dir);
    boostFS::create_directories(_proc_dir);
    boostFS::create_directories(_figure_dir);
}

