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

    EggnogDatabase(FileSystem* filesystem);
    ~EggnogDatabase();

    ERR_EGGNOG_DB download(EGGNOG_DB_TYPES type, std::string out_path);
    ERR_EGGNOG_DB open_sql(std::string& sql_path);
    std::string print_err();
    QuerySequence::EggnogResults get_eggnog_entry(std::string &accession);


private:
    const std::string FTP_EGGNOG_SQL  = "http://eggnogdb.embl.de/download/latest/eggnog-mapper-data/eggnog.db.gz";
    const std::string FTP_EGGNOG_DMND = "http://eggnogdb.embl.de/download/latest/eggnog-mapper-data/eggnog_proteins.dmnd.gz";
    const std::string FTP_EGGNOG_FASTA= "http://eggnogdb.embl.de/download/latest/eggnog-mapper-data/eggnog4.clustered_proteins.fa.gz";

    const std::string TEMP_SQL_GZ = "temp_egg_sql.gz";
    const std::string TEMP_DMND_GZ = "temp_egg_dmnd.gz";
    const std::string TEMP_FAST_GZ = "temp_egg_fasta.gz";

    // SQL Data
    const fp32        SQL_VERSION_CHANGE    = 1;            // Version in which SQL tables changed (just placeholder for now)
    const std::string SQL_MEMBER_TABLE      = "member";     // Changed to 'orthologs' in newer versions
    const std::string SQL_MEMBER_GROUP      = "groups";
    const std::string SQL_MEMBER_NAME       = "name";
    const std::string SQL_MEMBER_ORTHOINDEX = "orthoindex";
    const std::string SQL_MEMBER_PNAME      = "pname";
    const std::string SQL_MEMBER_GO         = "go";
    const std::string SQL_MEMBER_KEGG       = "kegg";
    const std::string SQL_EGGNOG_TABLE      = "eggnog";     // Newer versions use this table for annotation info
    const std::string SQL_EGGNOG_NAME       = "eggnog.name";
    const std::string SQL_EGGNOG_PNAME      = "seq.pname";
    const std::string SQL_EGGNOG_SEQ_NAME   = "seq.name";
    const std::string SQL_EGGNOG_GOS        = "gene_ontology.gos";
    const std::string SQL_EGGNOG_KEGG       = "kegg.ko";
    const std::string SQL_EGGNOG_BIGG       = "bigg.reaction";
    const std::string SQL_EVENT_TABLE       = "event";
    const std::string SQL_EVENT_LEVEL       = "level";
    const std::string SQL_EVENT_SIDE1       = "side1";
    const std::string SQL_EVENT_SIDE2       = "side2";
    const std::string SQL_EVENT_I           = "i";

    typename std::unordered_map<std::string, vect_str_t>::const_iterator _it_vect_str;
    SQLDatabaseHelper *_pSQLDatabase;
    FileSystem        *_pFilesystem;
    std::string        _err_msg;
    fp32               _version;
    static const std::unordered_map<std::string,std::string> EGGNOG_LEVELS;   // Mappings from tax lvl to full name
    static const std::unordered_map<std::string, vect_str_t> LEVEL_CONTENT;
    static const vect_str_t                                  TAXONOMIC_RESOLUTION;

    void get_tax_scope(std::string&, QuerySequence::EggnogResults&);
    void get_sql_data(QuerySequence::EggnogResults&, SQLDatabaseHelper&);
    std::string format_sql_data(std::string&);
    void get_og_query(QuerySequence::EggnogResults&);
    void get_member_ogs(QuerySequence::EggnogResults& eggnog_results);
    member_orthologs_t get_member_orthologs(member_orthologs_t &member_orthologs,
                              std::string &best_hit,
                              std::set<std::string> &target_lvls);
    void get_annotations(set_str_t& orthologs, QuerySequence::EggnogResults& eggnog_results);

};


#endif //ENTAP_EGGNOGDATABASE_H
