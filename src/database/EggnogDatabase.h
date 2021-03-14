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

#ifndef ENTAP_EGGNOGDATABASE_H
#define ENTAP_EGGNOGDATABASE_H

#include "../common.h"
#include "SQLDatabaseHelper.h"
#include "../FileSystem.h"
#include "../EntapGlobals.h"
#include "../QuerySequence.h"

class EggnogDatabase {

public:

    typedef std::map<std::string,set_str_t> member_orthologs_t;

    typedef enum {
        EGGNOG_SQL,         // EggNOG SQL database
        EGGNOG_DIAMOND,     // EggNOG DIAMOND database
        EGGNOG_FASTA        // EggNOG fasta (non-diamond)

    } EGGNOG_DB_TYPES;

    typedef enum {

        ERR_EGG_OK=0,
        ERR_EGG_SQL_FTP,
        ERR_EGG_SQL_DECOMP,
        ERR_EGG_SQL_FILE_EXISTS,
        ERR_EGG_SQL_OPEN,
        ERR_EGG_DMND_FTP,
        ERR_EGG_DMND_DECOMP,
        ERR_EGG_FASTA_FTP,
        ERR_EGG_FASTA_DECOMP,

    } ERR_EGGNOG_DB;

    typedef enum {
        EGGNOG_VERSION_UNKONWN=0,   // We could NOT find version and did not assume one
        EGGNOG_VERSION_EARLIER,     // EggNOG Database v4.1 (default)
        EGGNOG_VERSION_4_5_1,
        EGGNOG_VERSION_MAX
    } EGGNOG_SQL_VERSION;

    typedef enum {
        EGGNOG_DATA_GO=0,
        EGGNOG_DATA_KEGG,
        EGGNOG_DATA_BIGG,
        EGGNOG_DATA_PNAME
    } EGGNOG_DATA_TYPES;

    EggnogDatabase(FileSystem* filesystem, EntapDatabase* entap_data, QueryData* queryData);
    ~EggnogDatabase();

    ERR_EGGNOG_DB download(EGGNOG_DB_TYPES type, std::string out_path);
    ERR_EGGNOG_DB open_sql(std::string& sql_path);
    std::string print_err();
    void get_eggnog_entry(QuerySequence::EggnogResults *eg);


private:
    // EggNOG 4.1 URLs
    const std::string FTP_EGGNOG_SQL  = "http://eggnog5.embl.de/download/eggnog_4.1/eggnog-mapper-data/eggnog.db.gz";
    const std::string FTP_EGGNOG_DMND = "http://eggnog5.embl.de/download/eggnog_4.1/eggnog-mapper-data/eggnog_proteins.dmnd.gz";
    const std::string FTP_EGGNOG_FASTA= "http://eggnog5.embl.de/download/eggnog_4.1/eggnog-mapper-data/eggnog4.clustered_proteins.fa.gz";

    const std::string TEMP_SQL_GZ = "temp_egg_sql.gz";
    const std::string TEMP_DMND_GZ = "temp_egg_dmnd.gz";
    const std::string TEMP_FAST_GZ = "temp_egg_fasta.gz";

    // SQL Data Constants (from EggNOG SQL Database)

    /*      emapper.db-4.5.1
     * bigg
     *      name            VARCHAR(32)
     *      reaction        VARCHAR(32)
     * eggnog
     *      name            VARCHAR(32)
     *      group           TEXT
     * event
     *      i               INTEGER
     *      level           VARCHAR(16)
     *      og              VARCHAR(16)
     *      side1           TEXT
     *      side2           TEXT
     * gene_ontology
     *      name            VARCHAR(32)
     *      gos             TEXT
     * kegg
     *      name            VARCHAR(32)
     *      ko              VARCHAR(32)
     * og
     *      og              VARCHAR(16)
     *      level           VARCHAR(16)
     *      nm              INTEGER
     *      description     TEXT
     *      COG_categories  VARCHAR(8)
     *      GO_freq         TEXT
     *      KEGG_freq       TEXT
     *      SMART_freq      TEXT
     *      proteins        TEXT
     * orthologs
     *      name            VARCHAR(32)
     *      orthoindex      TEXT
     * seq
     *      name            VARCHAR(32)
     *      pname           VARCHAR(32)
     * version
     *      version         VARCHAR(16)
     *
     *
     */
    const std::string SQL_MEMBER_TABLE_1    = "orthologs";
    const std::string SQL_EGGNOG_TABLE      = "eggnog";     // Used for annotation info
    const std::string SQL_EGGNOG_KEGG       = "kegg.ko";
    const std::string SQL_EGGNOG_KEGG_NAME  = "kegg.name";
    const std::string SQL_EGGNOG_BIGG       = "bigg.reaction";
    const std::string SQL_EGGNOG_BIGG_NAME  = "bigg.name";
    const std::string SQL_EGGNOG_PNAME      = "seq.pname";
    const std::string SQL_EGGNOG_SEQ_NAME   = "seq.name";
    const std::string SQL_EGGNOG_GOS        = "gene_ontology.gos";
    const std::string SQL_EGGNOG_GO_NAME    = "gene_ontology.name";

    // emapper.db-earlier (currently being used)
    const std::string SQL_MEMBER_TABLE      = "member";     // Changed to 'orthologs' in newer versions

    // Both databases
    const std::string SQL_MEMBER_GROUP      = "groups";
    const std::string SQL_MEMBER_NAME       = "name";
    const std::string SQL_MEMBER_ORTHOINDEX = "orthoindex";
    const std::string SQL_MEMBER_PNAME      = "pname";
    const std::string SQL_MEMBER_GO         = "go";
    const std::string SQL_MEMBER_KEGG       = "kegg";
    const std::string SQL_EGGNOG_NAME       = "eggnog.name";

    const std::string SQL_EVENT_TABLE       = "event";
    const std::string SQL_EVENT_LEVEL       = "level";
    const std::string SQL_EVENT_SIDE1       = "side1";
    const std::string SQL_EVENT_SIDE2       = "side2";
    const std::string SQL_EVENT_I           = "i";

    typename std::unordered_map<std::string, vect_str_t>::const_iterator _it_vect_str;
    SQLDatabaseHelper *mpSQLDatabase;
    FileSystem        *mpFileSystem;
    EntapDatabase     *mpEntapDatabase;
    QueryData         *mpQueryData;         // Used to control header information
    std::string        mErrMsg;
    ERR_EGGNOG_DB      mErrCode;
    std::string        mSQLMemberTable;
    EGGNOG_SQL_VERSION mSQLVersion;
    uint16              mVersionMajor;
    uint16              mVersionMinor;
    uint16              mVersionRev;

    static const std::unordered_map<std::string,std::string> EGGNOG_LEVELS;   // Mappings from tax lvl to full name
    static const std::unordered_map<std::string, vect_str_t> LEVEL_CONTENT;
    static const vect_str_t                                  TAXONOMIC_RESOLUTION;

    void get_tax_scope(QuerySequence::EggnogResults*);
    void get_additional_sql_data(QuerySequence::EggnogResults* eggnogResults);
    std::string format_sql_data(std::string&);
    void get_og_query(QuerySequence::EggnogResults* eggnogResults);
    void get_member_ogs(QuerySequence::EggnogResults* eggnog_results);
    member_orthologs_t get_member_orthologs(member_orthologs_t &member_orthologs,
                              std::string &best_hit,
                              std::set<std::string> &target_lvls);
    void get_annotations(set_str_t& orthologs, QuerySequence::EggnogResults* eggnog_results);
    void set_error(std::string msg, ERR_EGGNOG_DB code);
    void set_database_version();
    void update_dataset(set_str_t &set, EGGNOG_DATA_TYPES datatype, std::string data);
};


#endif //ENTAP_EGGNOGDATABASE_H
