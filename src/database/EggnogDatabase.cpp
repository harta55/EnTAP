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

#include "EggnogDatabase.h"

const std::map<std::string,std::string> EGGNOG_LEVELS = {
        {"acoNOG", "Aconoidasida"},
        {"agaNOG", "Agaricales"},
        {"agarNOG", "Agaricomycetes"},
        {"meNOG", "Animals"},
        {"apiNOG", "Apicomplexa"},
        {"arthNOG", "Arthrodermataceae"},
        {"artNOG", "Arthropoda"},
        {"ascNOG", "Ascomycota"},
        {"aveNOG", "Aves"},
        {"basNOG", "Basidiomycota"},
        {"biNOG", "Bilateria"},
        {"braNOG", "Brassicales"},
        {"carNOG", "Carnivora"},
        {"chaNOG", "Chaetomiaceae"},
        {"chloroNOG", "Chlorophyta"},
        {"chorNOG", "Chordata"},
        {"chrNOG", "Chromadorea"},
        {"cocNOG", "Coccidia"},
        {"cryNOG", "Cryptosporidiidae"},
        {"debNOG", "Debaryomycetaceae"},
        {"dipNOG", "Diptera"},
        {"dotNOG", "Dothideomycetes"},
        {"droNOG", "Drosophilidae"},
        {"euNOG", "Eukaryotes"},
        {"eurotNOG", "Eurotiales"},
        {"euroNOG", "Eurotiomycetes"},
        {"fiNOG", "Fishes"},
        {"fuNOG", "Fungi"},
        {"haeNOG", "Haemosporida"},
        {"homNOG", "Hominidae"},
        {"hymNOG", "Hymenoptera"},
        {"hypNOG", "Hypocreales"},
        {"inNOG", "Insects"},
        {"kinNOG", "Kinetoplastida"},
        {"lepNOG", "Lepidoptera"},
        {"lilNOG", "Liliopsida"},
        {"magNOG", "Magnaporthales"},
        {"maNOG", "Mammals"},
        {"necNOG", "Nectriaceae"},
        {"nemNOG", "Nematodes"},
        {"onyNOG", "Onygenales"},
        {"opiNOG", "Opisthokonts"},
        {"perNOG", "Peronosporales"},
        {"pleNOG", "Pleosporales"},
        {"poaNOG", "Poales"},
        {"prNOG", "Primates"},
        {"rhaNOG", "Rhabditida"},
        {"roNOG", "Rodents"},
        {"sacNOG", "Saccharomycetaceae"},
        {"saccNOG", "Saccharomycetes"},
        {"sorNOG", "Sordariales"},
        {"sordNOG", "Sordariomycetes"},
        {"strNOG", "Streptophyta"},
        {"spriNOG", "Supraprimates"},
        {"treNOG", "Tremellales"},
        {"veNOG", "Vertebrates"},
        {"virNOG", "Viridiplantae"},
        {"bactNOG", "Bacteria"},
        {"bctoNOG", "Bacteroidetes"},
        {"chlNOG", "Chlorobi"},
        {"cyaNOG", "Cyanobacteria"},
        {"proNOG", "Proteobacteria"},
        {"gproNOG", "Proteobacteria_gamma"},
        {"firmNOG", "Firmicutes"},
        {"deiNOG", "Deinococcusthermus"},
        {"aproNOG", "Proteobacteria_alpha"},
        {"bproNOG", "Proteobacteria_beta"},
        {"dproNOG", "Proteobacteria_delta"},
        {"eproNOG", "Proteobacteria_epsilon"},
        {"chlorNOG", "Chloroflexi"},
        {"fusoNOG", "Fusobacteria"},
        {"aciNOG", "Acidobacteria"},
        {"delNOG", "delta/epsilon"},
        {"verNOG", "Verrucomicrobia"},
        {"bacNOG", "Bacilli"},
        {"flaNOG", "Flavobacteriia"},
        {"sphNOG", "Sphingobacteriia"},
        {"cloNOG", "Clostridia"},
        {"bacteNOG", "Bacteroidia"},
        {"aquNOG", "Aquificae"},
        {"chloNOG", "Chloroflexi"},
        {"therNOG", "Thermotogae"},
        {"defNOG", "Deferribacteres"},
        {"actNOG", "Actinobacteria"},
        {"verrNOG", "Verrucomicrobiae"},
        {"plaNOG", "Planctomycetes"},
        {"spiNOG", "Spirochaetes"},
        {"chlaNOG", "Chlamydiae"},
        {"acidNOG", "Acidobacteriia"},
        {"dehNOG", "Dehalococcoidetes"},
        {"synNOG", "Synergistetes"},
        {"eryNOG", "Erysipelotrichi"},
        {"tenNOG", "Tenericutes"},
        {"cytNOG", "Cytophagia"},
        {"negNOG", "Negativicutes"},
        {"arNOG", "Archaea"},
        {"creNOG", "Crenarchaeota"},
        {"eurNOG", "Euryarchaeota"},
        {"metNOG", "Methanobacteria"},
        {"methNOG", "Methanococci"},
        {"halNOG", "Halobacteria"},
        {"theNOG", "Thermoplasmata"},
        {"thermNOG", "Thermococci"},
        {"arcNOG", "Archaeoglobi"},
        {"methaNOG", "Methanomicrobia"},
        {"thaNOG", "Thaumarchaeota"},
        {"NOG", "Ancestor"}
};


EggnogDatabase::EggnogDatabase(FileSystem* filesystem) {
    _pFilesystem = filesystem;
    _pSQLDatabase = nullptr;
    _err_msg = "";
}

EggnogDatabase::~EggnogDatabase() {
    FS_dprint("Killing object - EggNOG Database");
    delete _pSQLDatabase;   // closes on SQLDatabaseHelper destructor
}

EggnogDatabase::ERR_EGGNOG_DB EggnogDatabase::download(EggnogDatabase::EGGNOG_DB_TYPES type, std::string out_path) {
    std::string temp_path;


    switch (type) {
        case EGGNOG_SQL:
            FS_dprint("Downloading EggNOG SQL database...");
            temp_path = PATHS(_pFilesystem->get_temp_outdir(), TEMP_SQL_GZ);
            if (!_pFilesystem->download_ftp_file(FTP_EGGNOG_SQL, temp_path)) {
                FS_dprint("Unable to download from FTP address at: " + FTP_EGGNOG_SQL);
                return ERR_EGG_SQL_FTP;
            }

            FS_dprint("Success! Decompressing...");

            if (!_pFilesystem->decompress_file(temp_path, out_path, FileSystem::FILE_GZ)) {
                FS_dprint("Unable to decompress file at: " + temp_path);
                return ERR_EGG_SQL_DECOMP;
            }
            FS_dprint("Success! EggNOG SQL database sent to: " + out_path);
            break;

        case EGGNOG_DIAMOND:
            FS_dprint("Downloading EggNOG DIAMOND database...");
            temp_path = PATHS(_pFilesystem->get_temp_outdir(), TEMP_DMND_GZ);
            if (!_pFilesystem->download_ftp_file(FTP_EGGNOG_DMND, temp_path)) {
                FS_dprint("Unable to download from FTP address at: " + FTP_EGGNOG_DMND);
                return ERR_EGG_DMND_FTP;
            }

            FS_dprint("Success! Decompressing...");

            if (!_pFilesystem->decompress_file(temp_path, out_path, FileSystem::FILE_GZ)) {
                FS_dprint("Unable to decompress file at: " + temp_path);
                _pFilesystem->delete_file(out_path);
                return ERR_EGG_DMND_DECOMP;
            }
            FS_dprint("Success! EggNOG DIAMOND database sent to: " + out_path);
            break;

        case EGGNOG_FASTA:
            FS_dprint("Downloading EggNOG DIAMOND database...");
            temp_path = PATHS(_pFilesystem->get_temp_outdir(), TEMP_FAST_GZ);
            if (!_pFilesystem->download_ftp_file(FTP_EGGNOG_FASTA, temp_path)) {
                FS_dprint("Unable to download from FTP address at: " + FTP_EGGNOG_FASTA);
                return ERR_EGG_FASTA_FTP;
            }

            FS_dprint("Success! Decompressing...");

            if (!_pFilesystem->decompress_file(temp_path, out_path, FileSystem::FILE_GZ)) {
                FS_dprint("Unable to decompress file at: " + temp_path);
                _pFilesystem->delete_file(out_path);
                return ERR_EGG_FASTA_DECOMP;
            }
            FS_dprint("Success! EggNOG FASTA sent to: " + out_path);
            break;

        default:
            break;

    }
    return ERR_EGG_OK;
}

EggnogDatabase::ERR_EGGNOG_DB EggnogDatabase::open_sql(std::string& sql_path) {
    FS_dprint("Opening SQL database...");
    if (!_pFilesystem->file_exists(sql_path)) {
        FS_dprint("File does not exist at: " + sql_path);
        return ERR_EGG_SQL_FILE_EXISTS;
    }

    if (_pSQLDatabase != nullptr) {
        FS_dprint("Database already created");
        return ERR_EGG_OK;
    }

    _pSQLDatabase = new SQLDatabaseHelper();
    if (!_pSQLDatabase->open(sql_path)) {
        FS_dprint("Unable to open SQL database");
        return ERR_EGG_SQL_OPEN;
    }

    return ERR_EGG_OK;
}

std::string EggnogDatabase::print_err() {
    return _err_msg;
}

QuerySequence::EggnogResults EggnogDatabase::get_eggnog_entry(std::string &accession) {
    QuerySequence::EggnogResults eggnog_data = {};

    eggnog_data.best_hit_query = accession;

    // Get member ogs (0A01R@biNOG,0V8CP@meNOG) from best hit query
    get_member_ogs(eggnog_data);



    return eggnog_data;
}



/**
 * ======================================================================
 * Function void EggnogDatabase::get_tax_scope(std::string &raw_scope,
                                    QuerySequence::EggnogResults &eggnogResults)
 *
 * Description          - Pulls the readable taxonomic scope from the EggnogLevels.h
 *                        file.
 *
 * Notes                - These may change, they are not pulled at configu (hardcoded)
 *
 * @param raw_scope     - Raw tax scope from EggNOG output (ie. virNOG[6])
 * @param eggnogResults - Current query sequence Eggnog struc
 *
 * @return              - None
 * ======================================================================
 */
void EggnogDatabase::get_tax_scope(std::string &raw_scope,
                              QuerySequence::EggnogResults &eggnogResults) {
    // Lookup/Assign Tax Scope

    if (!raw_scope.empty()) {
        uint16 p = (uint16) (raw_scope.find("NOG"));
        if (p != std::string::npos) {
            eggnogResults.tax_scope = raw_scope.substr(0,p+3);
            eggnogResults.tax_scope_readable = EGGNOG_LEVELS.at(eggnogResults.tax_scope);
            return;
        }
    }
    eggnogResults.tax_scope  = raw_scope;
    eggnogResults.tax_scope_readable = "";
}


/**
 * ======================================================================
 * Function void EggnogDatabase::get_sql_data(QuerySequence::EggnogResults &eggnogResults,
 *                                          DatabaseHelper &database)
 *
 * Description          - Query EggNOG SQL database and pull relevant info
 *                      - Sets values in EggnogResults struct
 *
 * Notes                - None
 *
 * @param database      - DatabaseHelper object of Eggnog database
 * @param eggnogResults - Current query sequence Eggnog struc
 *
 * @return              - None
 * ======================================================================
 */
void EggnogDatabase::get_sql_data(QuerySequence::EggnogResults &eggnogResults, SQLDatabaseHelper &database) {
    // Lookup description, KEGG, protein domain from SQL database
    if (!eggnogResults.og_key.empty()) {
        std::vector<std::vector<std::string>>results;
        std::string sql_kegg;
        std::string sql_desc;
        std::string sql_protein;

        char *query = sqlite3_mprintf(
                "SELECT description, KEGG_freq, SMART_freq FROM og WHERE og=%Q",
                eggnogResults.og_key.c_str());
        try {
            results = database.query(query);
            sql_desc = results[0][0];
            sql_kegg = results[0][1];
            sql_protein = results[0][2];
            if (!sql_desc.empty() && sql_desc.find("[]") != 0) eggnogResults.description = sql_desc;
            if (!sql_kegg.empty() && sql_kegg.find("[]") != 0) {
                eggnogResults.sql_kegg = format_sql_data(sql_kegg);
            }
            if (!sql_protein.empty() && sql_protein.find("{}") != 0){
                eggnogResults.protein_domains = format_sql_data(sql_protein);
            }
        } catch (std::exception &e) {
            // Do not fatal error
            FS_dprint(e.what());
        }
    }
}


/**
 * ======================================================================
 * Function void EggnogDatabase::get_og_query(QuerySequence::EggnogResults &eggnogResults)
 *
 * Description          - Find specific OG that aligned to database (based on tax scope)
 *                      - Sets values in EggnogResults struct
 *
 * Notes                - None
 *
 * @param eggnogResults - Current query sequence Eggnog struc
 *
 * @return              - None
 * ======================================================================
 */
void EggnogDatabase::get_og_query(QuerySequence::EggnogResults &eggnogResults) {
    // Find OG query was assigned to
    std::string temp;
    if (!eggnogResults.ogs.empty()) {
        std::istringstream ss(eggnogResults.ogs);
        std::unordered_map<std::string,std::string> og_map; // Not fully used right now
        while (std::getline(ss,temp,',')) {
            uint16 p = (uint16) temp.find("@");
            og_map[temp.substr(p+1)] = temp.substr(0,p);
        }
        eggnogResults.og_key = "";
        if (og_map.find(eggnogResults.tax_scope) != og_map.end()) {
            eggnogResults.og_key = og_map[eggnogResults.tax_scope];
        }
    }
}


/**
 * ======================================================================
 * Function std::string EggnogDatabase::format_sql_data(std::string &input)
 *
 * Description          - Parses SQL data and puts it in generic format for EnTAP
 *
 * Notes                - None
 *
 * @param input         - Unformatted SQL data
 *
 * @return              - Formatted string
 * ======================================================================
 */
std::string EggnogDatabase::format_sql_data(std::string &input) {
    enum FLAGS {
        DOM    = 0x01,
        INNER  = 0x02,
        INNER2 = 0x04,
        STR    = 0x08,
        FOUND  = 0x10
    };

    unsigned char bracketFlag = 0x0;
    std::string output = "";

    for (char c : input) {
        if (c == '{') {
            bracketFlag |= DOM;
        } else if (c == '[' && !(bracketFlag & INNER)) {
            bracketFlag |= INNER;
            if (bracketFlag & DOM) output += " (";
        } else if (c == '[' && !(bracketFlag & INNER2)) {
            bracketFlag |= INNER2;
        } else if (c == '}') {
            bracketFlag &= ~DOM;
            bracketFlag &= ~FOUND;
            output = output.substr(0,output.length()-2); // Remove trailing ', '
        } else if (c == ']' && (bracketFlag & INNER2)) {
            bracketFlag &= ~INNER2;
            bracketFlag &= ~FOUND;
            output += ", ";
        } else if (c == ']' && (bracketFlag & INNER)) {
            bracketFlag &= ~INNER;
            bracketFlag &= ~FOUND;
            output = output.substr(0,output.length()-2); // Remove trailing ', '
            if (bracketFlag & DOM) output += "), ";
        } else if (c == '\"') {
            bracketFlag ^= STR;
            if (!(bracketFlag & STR)) bracketFlag |= FOUND;
            if (!(bracketFlag & INNER)) bracketFlag &= ~FOUND;
        } else {
            if (bracketFlag & FOUND) continue;
            if (bracketFlag & STR) output += c;
        }
    }
    return output;
}


void EggnogDatabase::get_member_ogs(QuerySequence::EggnogResults& eggnog_results) {
    std::vector<std::vector<std::string>>results;

    if (eggnog_results.best_hit_query.empty()) return;

    char *query = sqlite3_mprintf(
            "SELECT %q FROM %q WHERE %q=%q",
            SQL_MEMBER_GROUP,
            SQL_MEMBER_TABLE,
            SQL_MEMBER_NAME,
            eggnog_results.best_hit_query);
    results = _pSQLDatabase->query(query);
    if (!results.empty()) {
        eggnog_results.member_ogs = results[0][1];
    }
}


