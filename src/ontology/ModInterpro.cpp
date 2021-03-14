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


//*********************** Includes *****************************
#include <csv.h>
#include "ModInterpro.h"
#include "../ExceptionHandler.h"
#include "../QueryData.h"


// Used for XML parsing
#if 0
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
using boost::property_tree::ptree;
#endif
#include "../FileSystem.h"
#include "../TerminalCommands.h"

//**************************************************************

std::vector<ENTAP_HEADERS> ModInterpro::DEFAULT_HEADERS = {
        ENTAP_HEADER_ONT_INTER_GO_BIO,
        ENTAP_HEADER_ONT_INTER_GO_CELL,
        ENTAP_HEADER_ONT_INTER_GO_MOLE,
        ENTAP_HEADER_ONT_INTER_PATHWAYS,
        ENTAP_HEADER_ONT_INTER_INTERPRO,
        ENTAP_HEADER_ONT_INTER_DATA_TYPE,
        ENTAP_HEADER_ONT_INTER_DATA_TERM,
        ENTAP_HEADER_ONT_INTER_EVAL
};

const std::vector<std::string> ModInterpro::INTERPRO_DATABASES ({
            "tigrfam",
            "sfld",
            "prodom",
            "hamap",
            "pfam",
            "smart",
            "cdd",
            "prositeprofiles",
            "prositepatterns",
            "superfamily",
            "prints",
            "panther",
            "gene3d",
            "pirsf",
            "coils",
            "mobidlite"
});

const std::string ModInterpro::INTERPRO_DEFAULT = "pfam";


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
EntapModule::ModVerifyData ModInterpro::verify_files() {
    ModVerifyData modVerifyData;
    std::string filename;

    filename        = INTERPRO_OUTPUT;
    mFinalBasepath = PATHS(mModOutDir, filename);
    filename += INTERPRO_EXT_TSV;
    mFinalOutpath = PATHS(mModOutDir, filename);
    modVerifyData.files_exist = mpFileSystem->file_exists(mFinalOutpath);
    return modVerifyData;
}

void ModInterpro::execute() {

    std::string interpro_cmd;
    std::string std_out;
    std::string blast;
    std::string temp_dir;
    int32       err_code;
    TerminalData      terminalData;

    mBlastp ? blast = PROTEIN_TAG : blast = NUCLEO_TAG;
    temp_dir = PATHS(mModOutDir, INTERPRO_TEMP);

    std_out      = PATHS(mModOutDir, INTERPRO_STD_OUT);
    interpro_cmd =
            mExePath    +
            " -i "          + mInputTranscriptome +
            " -b "          + mFinalBasepath +
            FLAG_SEQTYPE    + " " + blast     +
            FLAG_TEMP       + " " + temp_dir  +
            FLAG_GOTERM     +
            FLAG_IPRLOOK    +
            FLAG_PATHWAY;

    if (!mDatabases.empty()) {
        for (std::string &val : mDatabases) interpro_cmd += " --appl " + val;
    } else {
        throw ExceptionHandler("No InterPro databases selected!",
                ERR_ENTAP_RUN_INTERPRO);
    }

    // Execute command
    terminalData.command       = interpro_cmd;
    terminalData.print_files   = true;
    terminalData.suppress_std_err = false;
    terminalData.base_std_path = std_out;


    err_code = TC_execute_cmd(terminalData);
    if (err_code != 0) {
        mpFileSystem->delete_file(mFinalOutpath);
        throw ExceptionHandler("Error executing InterProScan\nInterProScan Error:\n"+ terminalData.err_stream,
                               ERR_ENTAP_RUN_INTERPRO);
    } else {
        mpFileSystem->delete_dir(temp_dir);
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
    std::string                           path_no_hits_faa;
    std::string                           path_no_hits_fnn;
    std::string                           path_hits_faa;
    std::string                           path_hits_fnn;
    std::map<std::string,InterProData>    interpro_map;
    go_format_t                           go_terms_parsed;
    uint32                                count_hits=0;
    uint32                                count_no_hits=0;

    FS_dprint("Beginning to parse InterProScan data...");
    if (mpFileSystem->file_exists(mFinalOutpath)) {
        FS_dprint("File found at: " + mFinalOutpath + " parsing...");
    } else {
        throw ExceptionHandler("Unable to locate InterProScan file at: " +
                               mFinalOutpath, ERR_ENTAP_PARSE_INTERPRO);
    }
    try {
        interpro_map = parse_tsv();
    } catch (ExceptionHandler const &e) {throw e;}

    FS_dprint("Success! Beginning to update query sequences...");

    path_no_hits_faa = PATHS(mProcDir, OUT_NO_HITS_FAA);
    path_no_hits_fnn = PATHS(mProcDir, OUT_NO_HITS_FNN);
    path_hits_faa    = PATHS(mProcDir, OUT_HITS_FAA);
    path_hits_fnn    = PATHS(mProcDir, OUT_HITS_FNN);

    std::ofstream file_no_hits_faa(path_no_hits_faa, std::ios::out | std::ios::app);
    std::ofstream file_no_hits_fnn(path_no_hits_fnn, std::ios::out | std::ios::app);
    std::ofstream file_hits_faa(path_hits_faa, std::ios::out | std::ios::app);
    std::ofstream file_hits_fnn(path_hits_fnn, std::ios::out | std::ios::app);

    if (!file_no_hits_faa.is_open() || !file_no_hits_fnn.is_open() ||
        !file_hits_faa.is_open()    || !file_hits_fnn.is_open()) {
        throw ExceptionHandler("Unable to open files for writing", ERR_ENTAP_FILE_IO);
    }

    // TODO stats
    QuerySequence::InterProResults interProResults;
    try {
        for (auto &pair : *mpQueryData->get_sequences_ptr()) {
            std::map<std::string, InterProData>::iterator it = interpro_map.find(pair.first);
            if (it != interpro_map.end()) {
                count_hits++;

                go_terms_parsed = mpEntapDatabase->format_go_delim(it->second.go_terms, ',');

                // Compile data TODO change...
                interProResults.e_value = float_to_sci(it->second.eval,2);
                interProResults.database_desc_id = it->second.databaseID + "(" + it->second.databaseDesc + ")";
                interProResults.database_type    = it->second.databasetype;
                interProResults.parsed_go        = go_terms_parsed;
                interProResults.interpro_desc_id = it->second.interID + "(" + it->second.interDesc + ")";
                interProResults.pathways         = it->second.pathways;
                interProResults.e_value_raw = it->second.eval;

                pair.second->add_alignment(mExecutionState, mSoftwareFlag, interProResults, INTERPRO_DATABASE_FLAG);

                if (!pair.second->get_sequence_n().empty()) file_hits_fnn << pair.second->get_sequence_n() << std::endl;
                if (!pair.second->get_sequence_p().empty()) file_hits_faa << pair.second->get_sequence_p() << std::endl;
            } else {
                // Not InterPro hit
                count_no_hits++;
                if (!pair.second->get_sequence_n().empty()) file_no_hits_fnn << pair.second->get_sequence_n() << std::endl;
                if (!pair.second->get_sequence_p().empty()) file_no_hits_faa << pair.second->get_sequence_p() << std::endl;
            }
        }
    } catch (std::exception &e) {
        throw ExceptionHandler(e.what(), ERR_ENTAP_PARSE_INTERPRO);
    }

    file_no_hits_faa.close();
    file_no_hits_fnn.close();
    file_hits_faa.close();
    file_hits_fnn.close();

    FS_dprint("Success! Calculating statistics...");
    mpFileSystem->format_stat_stream(stats_stream, "Ontology - InterProScan");
    stats_stream << "InterProScan statistics coming soon!";

    if (count_hits == 0) {
        stats_stream << "\nWarning: No InterProScan results found!";
    }
    stats_out = stats_stream.str();
    mpFileSystem->print_stats(stats_out);
    FS_dprint("Success! InterProScan finished");
}


#if 0

/**
* ======================================================================
* Function std::map<std::string,ModInterpro::InterProData> ModInterpro::parse_xml(void)
*
* Description          - Parses protein XML data from InterPro through boost library
*
* Notes                - Uses boost property tree
*
*
* @return              - Map of interpro data keyed to sequence id
*
* =====================================================================
*/
std::map<std::string,ModInterpro::InterProData> ModInterpro::parse_xml(void) {
    std::string                           seq_id;
    fp64                                  e_val;
    bool                                  inter;
    std::string                           pathway;
    ptree                                 pt;
    std::map<std::string,InterProData>    interpro_map;

    boost::property_tree::read_xml(mFinalOutpath, pt);
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
    return interpro_map;
}

#endif

/**
* ======================================================================
* Function std::map<std::string,ModInterpro::InterProData> ModInterpro::parse_tsv(void)
*
* Description          - Parses nucleotide data from interpro and returns
*                        parsed data map of BEST HITS according to E_value
*
* Notes                - None
*
* @return              - Map of parsed data keyed to sequence id
*
* =====================================================================
*/
std::map<std::string,ModInterpro::InterProData> ModInterpro::parse_tsv() {
    // Made them separate handles for xml vs tsv, can change though...might not need all info
    InterProData                        interProData;
    std::map<std::string,InterProData>  interpro_map;
    std::string temp_file_path;
    std::string query;
    std::string md5;
    std::string length;
    std::string database;
    std::string database_id;
    std::string database_desc;
    std::string start;
    std::string stop;
    fp64        eval;
    std::string status;
    std::string date;
    std::string interpro_id;
    std::string interpro_desc;
    std::string go_terms;       // GO:43112|GO:43111 format
    std::string pathways;       // KEGG: 00290+1.1.1.86|KEGG: 00770+1.1.1.86

    temp_file_path = format_interpro();

    io::CSVReader<INTERPRO_COL_NUM, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(temp_file_path);
    while (in.read_row(query, md5, length, database, database_id, database_desc,
                       start, stop, eval, status, date, interpro_id, interpro_desc,
                       go_terms, pathways)) {
        if (query.empty()) continue;

        std::map<std::string,InterProData>::iterator it = interpro_map.find(query);
        if (it != interpro_map.end() && it->second.eval < eval) continue;
        // Current hit is better
        if (!go_terms.empty() && go_terms.find('|') != std::string::npos) {
            std::replace(go_terms.begin(), go_terms.end(), '|', ',');
        }

        if (mpQueryData->get_sequence(query) == nullptr) {
            // InterProScan5 adds underscore to identifier
            if (query.find_last_of('_') != std::string::npos) {
                query = query.substr(0,query.find_last_of('_'));
                // Check if query is good
                if (mpQueryData->get_sequence(query) == nullptr) {
                    throw ExceptionHandler("InterPro Query: " + query + " not found in transcriptome!",
                                           ERR_ENTAP_PARSE_INTERPRO);
                }
            }
        }

        interProData = {};
        interProData.interID      = interpro_id;
        interProData.interDesc    = interpro_desc;
        interProData.databaseID   = database_id;
        interProData.databasetype = database;
        interProData.databaseDesc = database_desc;
        interProData.pathways     = pathways;
        interProData.go_terms     = go_terms;
        interProData.eval         = eval;
        interpro_map[query] = interProData;
    }

    mpFileSystem->delete_file(temp_file_path);
    return interpro_map;
}


/**
* ======================================================================
* Function std::string ModInterpro::format_interpro(void)
*
* Description          - Formats tsv file to be read easier with current
*                        lib (complains if tabs are not perfect)
*                      - This will add in tabs to keep everything consistent
*
*
* Notes                - None
*
* @return              - Path of temporary altered file
*
* =====================================================================
*/
std::string ModInterpro::format_interpro() {
    // Replace
    std::string path_temp;
    std::string line;
    uint16      tab_ct;

    path_temp = mFinalOutpath + "_temp";
    mpFileSystem->delete_file(path_temp);
    std::ifstream file_in(mFinalOutpath);
    std::ofstream file_temp(path_temp, std::ios::out | std::ios::app);
    while(std::getline(file_in, line)) {
        if (line.empty()) continue;
        file_temp << line;
        tab_ct = (uint16)std::count(line.begin(), line.end(), '\t');
        while (tab_ct < INTERPRO_COL_NUM - 1) {
            file_temp<<'\t';
            tab_ct++;
        }
        file_temp << std::endl;
    }
    file_in.close();
    file_temp.close();
    return path_temp;
}

ModInterpro::~ModInterpro() {
    FS_dprint("Killing object - ModInterpro");
}

bool ModInterpro::is_executable(std::string &exe) {
    std::string test_command;
    TerminalData terminalData;

    test_command = exe + " --version";

    terminalData.command = test_command;
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;

    return TC_execute_cmd(terminalData) == 0;
}

bool ModInterpro::valid_input(UserInput* userInput) {
    ent_input_multi_str_t databases;

    if (!userInput->has_input(INPUT_FLAG_INTERPRO)) return false;
    databases = userInput->get_user_input<ent_input_multi_str_t>(INPUT_FLAG_INTERPRO);

    for (std::string &data : databases) {
        LOWERCASE(data);
        if (!(FIND_VECT(data, INTERPRO_DATABASES))) {
            return false;
        }
    }
    return true;
}

ModInterpro::ModInterpro(std::string &ont, std::string &in, EntapDataPtrs& entap_data)
    : AbstractOntology(in, ont, entap_data, "InterProScan", DEFAULT_HEADERS){
    FS_dprint("Spawn Object - InterPro");

    mDatabases    = mpUserInput->get_user_input<ent_input_multi_str_t>(INPUT_FLAG_INTERPRO);
    mSoftwareFlag = ONT_INTERPRO_SCAN;
    mExePath      = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_INTERPRO_EXE);
}

bool ModInterpro::set_version() {
    return false;
}
