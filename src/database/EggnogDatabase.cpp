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

#include "EggnogDatabase.h"
#include "../QueryData.h"

const std::unordered_map<std::string,std::string> EggnogDatabase::EGGNOG_LEVELS = {
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

const std::unordered_map<std::string, vect_str_t> EggnogDatabase::LEVEL_CONTENT = {
        {"NOG", {"arNOG","bactNOG","euNOG","thaNOG","eurNOG","creNOG","synNOG","spiNOG","firmNOG","fusoNOG","aquNOG","aciNOG",
                "therNOG","tenNOG","proNOG","defNOG","plaNOG","actNOG","chloNOG","cyaNOG","deiNOG","bctoNOG","chlNOG",
                 "chlaNOG","verNOG","kinNOG","perNOG","virNOG","apiNOG","opiNOG","arcNOG","metNOG","methNOG","thermNOG",
                 "methaNOG","halNOG","theNOG","negNOG","cloNOG","eryNOG","bacNOG","acidNOG","delNOG","gproNOG","aproNOG",
                 "bproNOG","chlorNOG","dehNOG","cytNOG","bacteNOG","sphNOG","flaNOG","verrNOG","strNOG","chloroNOG",
                 "acoNOG","cocNOG","meNOG","fuNOG","dproNOG","eproNOG","braNOG","lilNOG","haeNOG","cryNOG","biNOG",
                 "basNOG","ascNOG","poaNOG","nemNOG","artNOG","chorNOG","agarNOG","treNOG","saccNOG","euroNOG","sordNOG",
                 "dotNOG","chrNOG","inNOG","veNOG","agaNOG","sacNOG","debNOG","eurotNOG","onyNOG","hypNOG","magNOG","sorNOG",
                 "pleNOG","rhaNOG","lepNOG","dipNOG","hymNOG","fiNOG","aveNOG","maNOG","arthNOG","necNOG","chaNOG","droNOG",
                 "spriNOG","carNOG","prNOG","roNOG","homNOG"}},
        {"arNOG", {"thaNOG","eurNOG","creNOG","arcNOG","metNOG","methNOG","thermNOG",
                   "methaNOG","halNOG","theNOG"}},
        {"bactNOG", {"synNOG","spiNOG","firmNOG","fusoNOG","aquNOG","aciNOG",
                     "therNOG","tenNOG","proNOG","defNOG","plaNOG","actNOG",
                     "chloNOG","cyaNOG","deiNOG","bctoNOG","chlNOG","chlaNOG",
                     "verNOG","negNOG","cloNOG","eryNOG","bacNOG","acidNOG",
                     "delNOG","gproNOG","aproNOG","bproNOG","chlorNOG","dehNOG",
                     "cytNOG","bacteNOG","sphNOG","flaNOG","verrNOG","dproNOG","eproNOG"}},
        {"euNOG", {"kinNOG","perNOG","virNOG","apiNOG","opiNOG","strNOG",
                   "chloroNOG","acoNOG","cocNOG","meNOG","fuNOG","braNOG",
                   "lilNOG","haeNOG","cryNOG","biNOG","basNOG","ascNOG",
                   "poaNOG","nemNOG","artNOG","chorNOG","agarNOG","treNOG",
                   "saccNOG","euroNOG","sordNOG","dotNOG","chrNOG","inNOG",
                   "veNOG","agaNOG","sacNOG","debNOG","eurotNOG","onyNOG",
                   "hypNOG","magNOG","sorNOG","pleNOG","rhaNOG","lepNOG",
                   "dipNOG","hymNOG","fiNOG","aveNOG","maNOG","arthNOG",
                   "necNOG","chaNOG","droNOG","spriNOG","carNOG","prNOG",
                   "roNOG","homNOG"}},
        {"eurNOG", {"arcNOG","metNOG","methNOG","thermNOG","methaNOG","halNOG","theNOG"}},
        {"firmNOG", {"negNOG","cloNOG","eryNOG","bacNOG"}},
        {"aciNOG", {"acidNOG"}},
        {"proNOG", {"delNOG","gproNOG","aproNOG","bproNOG","dproNOG","eproNOG"}},
        {"chloNOG", {"chlorNOG","dehNOG"}},
        {"bctoNOG", {"cytNOG","bacteNOG","sphNOG","flaNOG"}},
        {"verNOG", {"verrNOG"}},
        {"virNOG", {"strNOG","chloroNOG","braNOG","lilNOG","poaNOG"}},
        {"apiNOG", {"acoNOG","cocNOG","haeNOG","cryNOG"}},
        {"opiNOG", {"meNOG","fuNOG","biNOG","basNOG","ascNOG","nemNOG",
                    "artNOG","chorNOG","agarNOG","treNOG","saccNOG",
                    "euroNOG","sordNOG","dotNOG","chrNOG","inNOG",
                    "veNOG","agaNOG","sacNOG","debNOG","eurotNOG",
                    "onyNOG","hypNOG","magNOG","sorNOG","pleNOG","rhaNOG",
                    "lepNOG","dipNOG","hymNOG","fiNOG","aveNOG","maNOG",
                    "arthNOG","necNOG","chaNOG","droNOG","spriNOG","carNOG",
                    "prNOG","roNOG","homNOG"}},
        {"delNOG", {"dproNOG","eproNOG"}},
        {"strNOG", {"braNOG","lilNOG","poaNOG"}},
        {"acoNOG", {"haeNOG"}},
        {"cocNOG", {"cryNOG"}},
        {"meNOG", {"biNOG","nemNOG","artNOG","chorNOG","chrNOG","inNOG",
                   "veNOG","rhaNOG","lepNOG","dipNOG","hymNOG","fiNOG",
                   "aveNOG","maNOG","droNOG","spriNOG","carNOG","prNOG",
                   "roNOG","homNOG"}},
        {"fuNOG", {"basNOG","ascNOG","agarNOG","treNOG","saccNOG",
                   "euroNOG","sordNOG","dotNOG","agaNOG","sacNOG",
                   "debNOG","eurotNOG","onyNOG","hypNOG","magNOG",
                   "sorNOG","pleNOG","arthNOG","necNOG","chaNOG"}},
        {"lilNOG", {"poaNOG"}},
        {"biNOG", {"nemNOG","artNOG","chorNOG","chrNOG","inNOG",
                   "veNOG","rhaNOG","lepNOG","dipNOG","hymNOG",
                   "fiNOG","aveNOG","maNOG","droNOG","spriNOG",
                   "carNOG","prNOG","roNOG","homNOG"}},
        {"basNOG", {"agarNOG","treNOG","agaNOG"}},
        {"ascNOG", {"saccNOG","euroNOG","sordNOG",
                    "dotNOG","sacNOG","debNOG","eurotNOG","onyNOG",
                    "hypNOG","magNOG","sorNOG","pleNOG","arthNOG",
                    "necNOG","chaNOG"}},
        {"nemNOG", {"chrNOG","rhaNOG"}},
        {"artNOG", {"inNOG","lepNOG","dipNOG","hymNOG","droNOG"}},
        {"chorNOG", {"veNOG","fiNOG","aveNOG","maNOG","spriNOG",
                     "carNOG","prNOG","roNOG","homNOG"}},
        {"agarNOG", {"agaNOG"}},
        {"saccNOG", {"sacNOG","debNOG"}},
        {"euroNOG", {"eurotNOG","onyNOG","arthNOG"}},
        {"sordNOG", {"hypNOG","magNOG","sorNOG","necNOG","chaNOG"}},
        {"dotNOG", {"pleNOG"}},
        {"chrNOG", {"rhaNOG"}},
        {"inNOG", {"lepNOG","dipNOG","hymNOG","droNOG"}},
        {"veNOG", {"fiNOG","aveNOG","maNOG","spriNOG","carNOG","prNOG",
                   "roNOG","homNOG"}},
        {"onyNOG", {"arthNOG"}},
        {"hypNOG", {"necNOG"}},
        {"sorNOG", {"chaNOG"}},
        {"dipNOG", {"droNOG"}},
        {"maNOG", {"spriNOG","carNOG","prNOG","roNOG","homNOG"}},
        {"spriNOG", {"prNOG","roNOG","homNOG"}},
        {"prNOG", {"homNOG"}}
};

const vect_str_t EggnogDatabase::TAXONOMIC_RESOLUTION = {"apiNOG", "virNOG", "nemNOG", "artNOG",
                                         "maNOG","fiNOG", "aveNOG", "meNOG",
                                         "fuNOG", "opiNOG", "euNOG", "arNOG", "bactNOG",
                                         "NOG"};


EggnogDatabase::EggnogDatabase(FileSystem* filesystem, EntapDatabase* entap_data, QueryData* queryData) {
    mpFileSystem = filesystem;
    mpSQLDatabase = nullptr;
    mpEntapDatabase = entap_data;
    mpQueryData = queryData;
    mErrMsg = "";
    mErrCode = ERR_EGG_OK;
    mVersionMajor = 0;
    mVersionMinor = 0;
    mVersionRev   = 0;
    mSQLVersion = EGGNOG_VERSION_UNKONWN;
}

EggnogDatabase::~EggnogDatabase() {
    FS_dprint("Killing object - EggNOG Database");
    delete mpSQLDatabase;   // closes on SQLDatabaseHelper destructor
}

EggnogDatabase::ERR_EGGNOG_DB EggnogDatabase::download(EggnogDatabase::EGGNOG_DB_TYPES type, std::string out_path) {
    std::string temp_path;


    switch (type) {
        case EGGNOG_SQL:
            FS_dprint("Downloading EggNOG SQL database...");
            temp_path = PATHS(mpFileSystem->get_temp_outdir(), TEMP_SQL_GZ);
            if (!mpFileSystem->download_ftp_file(FTP_EGGNOG_SQL, temp_path)) {
                set_error("Unable to download from FTP address at: " + FTP_EGGNOG_SQL +
                          mpFileSystem->get_error(), ERR_EGG_SQL_FTP);
                return ERR_EGG_SQL_FTP;
            }

            FS_dprint("Success! Decompressing...");

            if (!mpFileSystem->decompress_file(temp_path, out_path, FileSystem::ENT_FILE_GZ)) {
                set_error("Unable to decompress file at: " + temp_path + mpFileSystem->get_error(),
                          ERR_EGG_SQL_DECOMP);
                return ERR_EGG_SQL_DECOMP;
            }
            FS_dprint("Success! EggNOG SQL database sent to: " + out_path);
            break;

        case EGGNOG_DIAMOND:
            FS_dprint("Downloading EggNOG DIAMOND database...");
            temp_path = PATHS(mpFileSystem->get_temp_outdir(), TEMP_DMND_GZ);
            if (!mpFileSystem->download_ftp_file(FTP_EGGNOG_DMND, temp_path)) {
                set_error("Unable to download from FTP address at: " + FTP_EGGNOG_DMND + mpFileSystem->get_error(),
                        ERR_EGG_DMND_FTP);
                mpFileSystem->delete_file(temp_path);
                return ERR_EGG_DMND_FTP;
            }

            FS_dprint("Success! Decompressing...");

            if (!mpFileSystem->decompress_file(temp_path, out_path, FileSystem::ENT_FILE_GZ)) {
                set_error("Unable to decompress file at: " + temp_path +mpFileSystem->get_error(),
                        ERR_EGG_DMND_DECOMP);
                mpFileSystem->delete_file(out_path);
                return ERR_EGG_DMND_DECOMP;
            }
            FS_dprint("Success! EggNOG DIAMOND database sent to: " + out_path);
            break;

        case EGGNOG_FASTA:
            FS_dprint("Downloading EggNOG DIAMOND database...");
            temp_path = PATHS(mpFileSystem->get_temp_outdir(), TEMP_FAST_GZ);
            if (!mpFileSystem->download_ftp_file(FTP_EGGNOG_FASTA, temp_path)) {
                set_error("Unable to download from FTP address at: " + FTP_EGGNOG_FASTA + mpFileSystem->get_error(),
                        ERR_EGG_FASTA_FTP);
                mpFileSystem->delete_file(temp_path);
                return ERR_EGG_FASTA_FTP;
            }

            FS_dprint("Success! Decompressing...");

            if (!mpFileSystem->decompress_file(temp_path, out_path, FileSystem::ENT_FILE_GZ)) {
                set_error("Unable to decompress file at: " + temp_path + mpFileSystem->get_error(),
                        ERR_EGG_FASTA_DECOMP);
                mpFileSystem->delete_file(out_path);
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
    if (!mpFileSystem->file_exists(sql_path)) {
        FS_dprint("File does not exist at: " + sql_path);
        return ERR_EGG_SQL_FILE_EXISTS;
    }

    if (mpSQLDatabase != nullptr) {
        FS_dprint("Database already created");
        return ERR_EGG_OK;
    }

    mpSQLDatabase = new SQLDatabaseHelper();
    if (!mpSQLDatabase->open(sql_path)) {
        FS_dprint("Unable to open SQL database");
        return ERR_EGG_SQL_OPEN;
    }

    set_database_version();
    return ERR_EGG_OK;
}

std::string EggnogDatabase::print_err() {
    return "\nEggNOG Database Error: " + mErrMsg;
}

void EggnogDatabase::get_eggnog_entry(QuerySequence::EggnogResults *eggnog_data) {
    std::set<std::string> unique_groups;    // Unique member orthologous groups
    std::string           temp;
    member_orthologs_t    member_orthologs;
    std::set<std::string> level_set;
    set_str_t             orthologs;        // Selected from member orthologs

    // Get member orthologous groups (0A01R@biNOG,0V8CP@meNOG) from best hit query
    get_member_ogs(eggnog_data);
    if (eggnog_data->member_ogs.empty()) return;

    // get SQL data and query info
//    get_og_query(eggnog_data);      // populated og_key to be used as index into SQL
//    get_sql_data(eggnog_data);


    // Get unique tax groups (split "0V8CP@meNOG" to meNOG) and max level
    std::istringstream iss(eggnog_data->member_ogs);
    while(std::getline(iss, temp, ',')) {
        unique_groups.insert(temp.substr(temp.find("@")+1));    // add meNOG
    }
    // For default taxonomic scope (may want to allow user to change later)
    for (const std::string &level : EggnogDatabase::TAXONOMIC_RESOLUTION) {
        if (unique_groups.find(level) != unique_groups.end()) {
            _it_vect_str = LEVEL_CONTENT.find(level);
            if (_it_vect_str != LEVEL_CONTENT.end()) {
                std::copy(_it_vect_str->second.begin(),
                      _it_vect_str->second.end(),
                      std::inserter(level_set,level_set.end()));
            }
            level_set.insert(level);
            eggnog_data->tax_scope_lvl_max = level + "[" + std::to_string(level_set.size()) + "]";
            // Get tax scope readable
            get_tax_scope(eggnog_data);
            break;
        }
    }

    // Get all member orthologs
    member_orthologs = get_member_orthologs(member_orthologs, eggnog_data->seed_ortholog, level_set);
    orthologs = member_orthologs["all"];        // default, can change

    if (!orthologs.empty()) {
        get_annotations(orthologs, eggnog_data);    // Pull final annotations from database
        try {
            get_additional_sql_data(eggnog_data);       // Pull any additional info from database
        } catch (...) {
            ; // Do not fatal error for this, it is not supported on all databases and info may not
              // be needed by the user
        }
    }
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
void EggnogDatabase::get_tax_scope(QuerySequence::EggnogResults *eggnogResults) {
    // Lookup/Assign Tax Scope

    eggnogResults->tax_scope_readable = "";
    eggnogResults->tax_scope  = eggnogResults->tax_scope_lvl_max;

    if (!eggnogResults->tax_scope_lvl_max.empty()) {
        uint16 p = (uint16) (eggnogResults->tax_scope_lvl_max.find("NOG"));
        if (p != std::string::npos) {
            eggnogResults->tax_scope = eggnogResults->tax_scope_lvl_max.substr(0,p+3);
            eggnogResults->tax_scope_readable = EGGNOG_LEVELS.at(eggnogResults->tax_scope);
            return;
        }
    }
}


/**
 * ======================================================================
 * Function void EggnogDatabase::get_additional_sql_data(QuerySequence::EggnogResults &eggnogResults,
 *                                          DatabaseHelper &database)
 *
 * Description          - Query EggNOG SQL database and pull additional information
 *
 * Notes                - None
 *
 * @param database      - DatabaseHelper object of Eggnog database
 * @param eggnogResults - Current query sequence Eggnog struc
 *
 * @return              - None
 * ======================================================================
 */
void EggnogDatabase::get_additional_sql_data(QuerySequence::EggnogResults *eggnogResults) {
    // Lookup description, KEGG, protein domain from SQL database

    if (mSQLVersion != EGGNOG_VERSION_EARLIER) return;  // Only supported for earlier versions of SQL database currently

    get_og_query(eggnogResults);    // Will lookup og_key

    if (!eggnogResults->og_key.empty()) {
        std::vector<std::vector<std::string>>results;
        std::string sql_kegg;
        std::string sql_desc;
        std::string sql_protein;

        char *query = sqlite3_mprintf(
                "SELECT description, KEGG_freq, SMART_freq FROM og WHERE og=%Q",
                eggnogResults->og_key.c_str());
        try {
            results = mpSQLDatabase->query(query);
            sql_desc = results[0][0];
            sql_kegg = results[0][1];
            sql_protein = results[0][2];
            if (!sql_desc.empty() && sql_desc.find("[]") != 0) eggnogResults->description = sql_desc;
#if 0
            if (!sql_kegg.empty() && sql_kegg.find("[]") != 0) {
                eggnogResults->kegg = format_sql_data(sql_kegg);
            }
#endif
            if (!sql_protein.empty() && sql_protein.find("{}") != 0){
                eggnogResults->protein_domains = format_sql_data(sql_protein);
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
void EggnogDatabase::get_og_query(QuerySequence::EggnogResults *eggnogResults) {
    // Find OG query was assigned to
    std::string temp;
    if (!eggnogResults->member_ogs.empty()) {
        std::istringstream ss(eggnogResults->member_ogs);
        std::unordered_map<std::string,std::string> og_map; // Not fully used right now
        while (std::getline(ss,temp,',')) {
            uint16 p = (uint16) temp.find("@");
            og_map[temp.substr(p+1)] = temp.substr(0,p);
        }
        eggnogResults->og_key = "";
        if (og_map.find(eggnogResults->tax_scope) != og_map.end()) {
            eggnogResults->og_key = og_map[eggnogResults->tax_scope];
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

void EggnogDatabase::get_member_ogs(QuerySequence::EggnogResults *eggnog_results) {
    std::vector<std::vector<std::string>>results;
    char *query;

    if (eggnog_results->seed_ortholog.empty()) return;

    if (mSQLVersion == EGGNOG_VERSION_4_5_1) {
        // emapper.db-4.5.1
        query = sqlite3_mprintf(
                "SELECT %q FROM %q WHERE %q=%Q",
                SQL_MEMBER_GROUP.c_str(),
                SQL_EGGNOG_TABLE.c_str(),
                SQL_MEMBER_NAME.c_str(),
                eggnog_results->seed_ortholog.c_str()
        );
    } else {
        // Older versions
        query = sqlite3_mprintf(
                "SELECT %q FROM %q WHERE %q=%Q",
                SQL_MEMBER_GROUP.c_str(),
                mSQLMemberTable.c_str(),
                SQL_MEMBER_NAME.c_str(),
                eggnog_results->seed_ortholog.c_str());
    }

    results = mpSQLDatabase->query(query);
    if (!results.empty()) {
        eggnog_results->member_ogs = results[0][0];
    }
}

EggnogDatabase::member_orthologs_t EggnogDatabase::get_member_orthologs(EggnogDatabase::member_orthologs_t &member_orthologs,
                                          std::string &best_hit,
                                          std::set<std::string> &target_lvls) {
    std::string                     query_taxon;  // Tax number from match
    std::string                     event_indexes;
    std::set<std::string>           target_members;
    char*                           sql_query;
    SQLDatabaseHelper::query_struct sql_results;

    query_taxon = best_hit.substr(0, best_hit.find_first_of('.'));    // "34740"
    target_members.insert(best_hit);                                  // 34740.HMEL017225-PA

    sql_query = sqlite3_mprintf(
            "SELECT %q FROM %q WHERE %q=%Q",
            SQL_MEMBER_ORTHOINDEX.c_str(),
            mSQLMemberTable.c_str(),
            SQL_MEMBER_NAME.c_str(),
            best_hit.c_str());
    sql_results = mpSQLDatabase->query(sql_query);
    if (!sql_results.empty()) {
        event_indexes = sql_results[0][0];
    } else return member_orthologs_t();

    if (event_indexes.empty()) return member_orthologs_t();
    // Can specify levels as well here
    sql_query = sqlite3_mprintf(
            "SELECT %q, %q, %q FROM %q WHERE %q IN %s AND %q IN %s",
            SQL_EVENT_LEVEL.c_str(),
            SQL_EVENT_SIDE1.c_str(),
            SQL_EVENT_SIDE2.c_str(),
            SQL_EVENT_TABLE.c_str(),
            SQL_EVENT_I.c_str(),
            mpSQLDatabase->format_string(event_indexes,',').c_str(),
            SQL_EVENT_LEVEL.c_str(),
            mpSQLDatabase->format_container(target_lvls).c_str()
    );
    sql_results = mpSQLDatabase->query(sql_query);

    std::map<std::pair<std::string,set_str_t>,
            std::set<std::pair<std::string,set_str_t>>> ortholog_map;

    for (std::vector<std::string> &hit : sql_results) {
//        std::string* level = &hit[0];
        std::string* side1 = &hit[1];
        std::string* side2 = &hit[2];

        // Vector of tax, id pairs
        std::vector<pair_str_t> side1_pairs;     // '6238' , 'CBG18195']
        std::vector<pair_str_t> side2_pairs;

        // Convert string of hits to tax, id pair for side1
        for (std::string &temp : split_string(*side1, ',')) {
            uint16 index = (uint16)temp.find_first_of('.');
            side1_pairs.push_back(std::make_pair(
                    temp.substr(0, index),
                    temp.substr(index+1)
            ));
        }
        // Convert string of hits to tax, id pair for side2
        for (std::string &temp : split_string(*side2, ',')) {
            uint16 index = (uint16)temp.find_first_of('.');
            side1_pairs.push_back(std::make_pair(
                    temp.substr(0, index),
                    temp.substr(index+1)
            ));
        }

        // Can add check for target taxonomy here
        vect_str_t target_taxa;
        vect_str_t targets;
        std::string mid;
        std::unordered_map<std::string,set_str_t> by_sp1;
        std::unordered_map<std::string,set_str_t> by_sp2;

        std::map<std::pair<std::string,set_str_t>,
                std::set<std::pair<std::string,set_str_t>>>::iterator iter;

        for (auto &pair : side1_pairs) {
            if ( (target_taxa.empty()) ||
                 (std::find(target_taxa.begin(),target_taxa.end(), pair.first) != target_taxa.end()) ||
                 (query_taxon.compare(pair.first) == 0)) {
                mid = pair.first + "." + pair.second;
                by_sp1[pair.first].insert(mid);// Default of empty set here
            }
        }
        for (auto &pair : side2_pairs) {
            if ( (target_taxa.empty()) ||
                 (std::find(target_taxa.begin(),target_taxa.end(), pair.first) != target_taxa.end()) ||
                 (query_taxon.compare(pair.first) == 0)) {
                mid = pair.first + "." + pair.second;
                by_sp2[pair.first].insert(mid); // Default of empty set here
            }
        }

        // merge side1 coorthologs
        if (!target_taxa.empty()) {
            targets = target_taxa;
        } else if (!by_sp2.empty()) {
            for (std::unordered_map<std::string,set_str_t>::iterator it = by_sp2.begin();
                    it != by_sp2.end(); it++) {
                targets.push_back(it->first);
            }
        }
        for (auto &pair : by_sp1) {
            if (!target_members.empty() && !pair.second.empty()) {
                std::pair<std::string,set_str_t> key1 = std::make_pair(
                        pair.first, pair.second
                ); // sort pair.second
                std::pair<std::string,set_str_t> key2;
                for (std::string &sp2 : targets) {
                    if (by_sp2.find(sp2) == by_sp2.end()) continue;
                    set_str_t co2 = by_sp2[sp2];
                    key2 = std::make_pair(sp2, co2);
                }
                iter = ortholog_map.find(key1);
                if (iter != ortholog_map.end()) {
                    iter->second.insert(key2);       // default empty set
                } else {
                    ortholog_map.emplace(key1, std::set<std::pair<std::string,set_str_t>>());
                    ortholog_map.at(key1).insert(key2);
                }
            }
        }

        // merge side2 coorthologs
        if (!target_taxa.empty()) {
            targets = target_taxa;
        } else if (!by_sp1.empty()) {
            for (std::unordered_map<std::string,set_str_t>::iterator it = by_sp1.begin();
                 it != by_sp1.end(); it++) {
                targets.push_back(it->first);
            }
        }
        for (auto &pair : by_sp2) {
            if (!target_members.empty() && !pair.second.empty()) {
                std::pair<std::string,set_str_t> key1 = std::make_pair(
                        pair.first, pair.second
                ); // sort pair.second
                std::pair<std::string,set_str_t> key2;
                for (std::string &sp2 : targets) {
                    if (by_sp1.find(sp2) == by_sp1.end()) continue;
                    set_str_t co2 = by_sp1[sp2];
                    key2 = std::make_pair(sp2, co2);
                }
                iter = ortholog_map.find(key1);
                if (iter != ortholog_map.end()) {
                    iter->second.insert(key2);       // default empty set
                } else {
                    ortholog_map.emplace(key1, std::set<std::pair<std::string,set_str_t>>());
                    ortholog_map.at(key1).insert(key2);
                }
            }
        }
    }

    member_orthologs_t all_orthologs {
            {"one2one", set_str_t()},
            {"one2many", set_str_t()},
            {"many2many", set_str_t()},
            {"many2one", set_str_t()},
            {"all", set_str_t()}
    };

    std::string otype_prefix;
    std::string otype;
    for (auto &pair : ortholog_map) {
        if (pair.second.size() == 1) {
            otype_prefix = "one2";
        } else {
            otype_prefix = "many2";
        }
        all_orthologs["all"].insert(pair.first.second.begin(), pair.first.second.end());

        for (auto &pair2 : pair.second) {
            if (pair2.second.size() == 1) {
                otype = otype_prefix + "one";
            } else {
                otype = otype_prefix + "many";
            }
            all_orthologs[otype].insert(pair.first.second.begin(), pair.first.second.end());
            all_orthologs[otype].insert(pair2.second.begin(), pair2.second.end());
            all_orthologs["all"].insert(pair2.second.begin(), pair2.second.end());
        }
    }
    return all_orthologs;
}

void EggnogDatabase::get_annotations(set_str_t& orthologs, QuerySequence::EggnogResults* eggnog_results) {

    char*               sql_query;
    set_str_t           all_gos;
    set_str_t           all_kegg;
    set_str_t           all_pnames;
    Compair<std::string>             pname_counter;
    set_str_t           all_bigg;
    std::string         delim_list;
    SQLDatabaseHelper::query_struct sql_results;

    // This is different depending on version on eggnog using
    if (mSQLVersion == EGGNOG_VERSION_4_5_1) {
        sql_query = sqlite3_mprintf(
                "SELECT %q, %q, %q, %q, %q FROM %q "\
                "LEFT JOIN seq on %q = %q "\
                "LEFT JOIN gene_ontology on %q = %q "\
                "LEFT JOIN kegg on %q = %q "\
                "LEFT JOIN bigg on %q = %q "\
                "WHERE %q in %s",
                SQL_EGGNOG_NAME.c_str(),
                SQL_EGGNOG_PNAME.c_str(),
                SQL_EGGNOG_GOS.c_str(),
                SQL_EGGNOG_KEGG.c_str(),
                SQL_EGGNOG_BIGG.c_str(),
                SQL_EGGNOG_TABLE.c_str(),
                SQL_EGGNOG_SEQ_NAME.c_str(), SQL_EGGNOG_NAME.c_str(),
                SQL_EGGNOG_GO_NAME.c_str(), SQL_EGGNOG_NAME.c_str(),
                SQL_EGGNOG_KEGG_NAME.c_str(), SQL_EGGNOG_NAME.c_str(),
                SQL_EGGNOG_BIGG_NAME.c_str(), SQL_EGGNOG_NAME.c_str(),
                SQL_EGGNOG_NAME.c_str(), mpSQLDatabase->format_container(orthologs).c_str()
        );
    } else {
        // Older versions
        sql_query = sqlite3_mprintf(
                "SELECT %q, %q, %q, %q FROM %q WHERE %q IN %s",
                SQL_MEMBER_NAME.c_str(),
                SQL_MEMBER_PNAME.c_str(),
                SQL_MEMBER_GO.c_str(),
                SQL_MEMBER_KEGG.c_str(),
                mSQLMemberTable.c_str(),
                SQL_MEMBER_NAME.c_str(), mpSQLDatabase->format_container(orthologs).c_str()
        );
    }

    sql_results = mpSQLDatabase->query(sql_query);
    if (!sql_results.empty()) {
        for (vect_str_t &data : sql_results) {
            update_dataset(all_pnames, EGGNOG_DATA_PNAME, data[1]);
            pname_counter.add_value(data[1]);
            update_dataset(all_gos, EGGNOG_DATA_GO, data[2]);
            update_dataset(all_kegg, EGGNOG_DATA_KEGG, data[3]);
            if (mSQLVersion == EGGNOG_VERSION_4_5_1)
                update_dataset(all_bigg, EGGNOG_DATA_BIGG, data[4]);
        }
        eggnog_results->pname  = container_to_string<std::string>(all_pnames,",");
        if (!pname_counter.empty()) {
            pname_counter.sort(true);
            if (pname_counter._sorted[0].second >= 2) {
                eggnog_results->predicted_gene = pname_counter._sorted[0].first;
            }
        }

        delim_list = container_to_string<std::string>(all_gos, ",");
        eggnog_results->parsed_go = mpEntapDatabase->format_go_delim(delim_list,',');
        eggnog_results->kegg = container_to_string<std::string>(all_kegg, ",");
        if (mSQLVersion == EGGNOG_VERSION_4_5_1)
            eggnog_results->bigg = container_to_string<std::string>(all_bigg, ",");
    } else {
        eggnog_results->pname = "";
        eggnog_results->parsed_go = go_format_t();
        eggnog_results->kegg = "";
        eggnog_results->bigg = "";
    }
}

void EggnogDatabase::set_error(std::string msg, ERR_EGGNOG_DB code) {
    FS_dprint(msg);
    mErrMsg = msg;
    mErrCode = code;
}

void EggnogDatabase::set_database_version() {
    char *query;
    std::string version_str;
    uint8 indexl;
    uint8 indexr;

    FS_dprint("Getting EggNOG database SQL version...");

    if (mpSQLDatabase == nullptr) {
        FS_dprint("Could not find EggNOG version, NULL");
        mSQLVersion = EGGNOG_VERSION_UNKONWN;
        return;
    }

    query = sqlite3_mprintf(
            "SELECT version FROM version"
    );

    try {
        FS_dprint("Executing sql cmd: " + std::string(query));
        version_str = mpSQLDatabase->query(query)[0][0];
        if (version_str.empty())
            throw std::runtime_error("NULL version from query");

        indexl = (uint8) version_str.find('.');
        indexr = (uint8) version_str.rfind('.');

        mVersionMajor = (uint16) std::stoi(version_str.substr(0,indexl));
        mVersionMinor = (uint16) std::stoi(version_str.substr(indexl+1, (uint32) (indexr - indexl - 1)));
        mVersionRev   = (uint16) std::stoi(version_str.substr(indexr+1));

        if (mVersionMajor == 4 && mVersionMinor == 5 && mVersionRev == 1) {
            mSQLVersion = EGGNOG_VERSION_4_5_1;
            FS_dprint("SQL Version set to 4.5.1");

        } else if (mVersionMajor <= 4 ){
            FS_dprint("WARNING: SQL Version less than 4, setting to 'earlier'");
            mSQLVersion = EGGNOG_VERSION_EARLIER;
            if (mpQueryData != nullptr)
                mpQueryData->header_set(ENTAP_HEADER_ONT_EGG_BIGG, false);   // Not supported for older version
        }
        FS_dprint("Success!");

    } catch (...) {
        FS_dprint("WARNING: couldn't get SQL version.Setting to earlier version by default");
        mSQLVersion = EGGNOG_VERSION_EARLIER;
        if (mpQueryData != nullptr)
            mpQueryData->header_set(ENTAP_HEADER_ONT_EGG_BIGG, false);       // Not supported for older version
    }

    if (mSQLVersion == EGGNOG_VERSION_4_5_1) {
        mSQLMemberTable = SQL_MEMBER_TABLE_1;
    } else {
        mSQLMemberTable = SQL_MEMBER_TABLE;
    }
}

void EggnogDatabase::update_dataset(set_str_t &set, EGGNOG_DATA_TYPES datatype, std::string data) {
    if (data.empty()) return;

    uint16 index_start;
    uint16 index_end;
    std::string term;

    switch (datatype) {
        case EGGNOG_DATA_GO:
            // F|GO:0000166|IEA,F|GO:0001882|IEA,F|GO:0001883|IEA
            while (data.find("GO:") != std::string::npos) {
                index_start = (uint16) data.find("GO:");
                index_end = (uint16) data.find("|", index_start);
                term = data.substr(index_start,data.find("|", index_start) - index_start);
                set.insert(term);
                data.erase(0,index_end);
            }
            break;

        case EGGNOG_DATA_KEGG:
            set.insert(data);
            break;

        case EGGNOG_DATA_BIGG:
        {
            std::istringstream iss(data);
            while(std::getline(iss, term, ',')) {
                set.insert(term);
            }
            break;
        }


        case EGGNOG_DATA_PNAME:
            set.insert(data);
            break;

        default:
            return;
    }
}