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

#ifndef ENTAP_ENTAPDATABASE_H
#define ENTAP_ENTAPDATABASE_H

#include "../EntapGlobals.h"
#include "SQLDatabaseHelper.h"
#include "../FileSystem.h"

#ifdef USE_BOOST    // Include boost serialization headers
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#else              // Include cereal serialization libraries
#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/set.hpp>

#endif

struct  GoEntry {
    std::string go_id;
    std::string level;
    std::string category;
    std::string term;
    int16       level_int;

    static const std::string UNKNOWN_LVL_STR;
    static const int16 UNKNOWN_LVL;

#ifdef USE_BOOST
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const uint32 v) {
        ar&go_id;
        ar&level;
        ar&category;
        ar&term;
    }
#else
    // Use CEREAL for serialization (ha..puns)
    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(
                go_id, level_int, level, category, term);
    }
#endif

    GoEntry() {
        // Set as defaults if we cannot find information for it
        go_id = "";
        level = UNKNOWN_LVL_STR;
        level_int = UNKNOWN_LVL;
        category = UNKNOWN_LVL_STR;
        term = "";
    }

    bool is_empty() {
        return (term.empty() || term == UNKNOWN_LVL_STR);
    }

    // Don't really care about order, just don't want duplicates
    //  std::set helper
    bool operator < (const GoEntry& rhs) const {
        return this->go_id < rhs.go_id;
    }

    bool operator == (const GoEntry& rhs) const {
        return this->go_id == rhs.go_id;
    }
};

namespace std {
    template<>
    struct hash<GoEntry> {
        std::size_t operator()(const GoEntry &x) const {
            return std::hash<std::string>()(x.go_id);
        }
    };
}


struct TaxEntry {
    std::string tax_id;
    std::string lineage;
    std::string tax_name;

#ifdef USE_BOOST
    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive & ar, const uint32 v) {
        ar&tax_id;
        ar&lineage;
        ar&tax_name;
    }
#else
    // Use CEREAL for serialization
    template<class Archive>
    void serialize(Archive & archive) {
        archive(
                tax_id, lineage, tax_name);
    }
#endif

    bool is_empty() {
        return this->tax_id.empty() && this->lineage.empty();
    }
    TaxEntry() {
        tax_id = "";
        lineage = "";
        tax_name = "";
    }
};

typedef std::set<GoEntry> go_format_t;

struct UniprotEntry {
    std::string database_x_refs;        // DR   OrthoDB; VOG090000I8; -.
    std::string comments;               // CC   -!- FUNCTION: Transcription activation. {ECO:0000305}.
    std::string uniprot_id;             // ID   001R_FRG3G              Reviewed;         256 AA.
    go_format_t go_terms;
    std::string kegg_terms;
#ifdef USE_BOOST
    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive & ar, const uint32 v) {
        ar&uniprot_id;
        ar&database_x_refs;
        ar&comments;
        ar&go_terms;
        ar&kegg_terms;
    }
#else
    // Use CEREAL for serialization
    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(
                database_x_refs, comments, uniprot_id, go_terms, kegg_terms);
    }
#endif

    UniprotEntry() {
        database_x_refs = "";
        comments = "";
        uniprot_id = "";
        kegg_terms = "";
        go_terms = {};
    }

    std::string print(void) {
        std::stringstream out;

        out <<
            "UniProt ID: " << uniprot_id <<
            "XRefs: "      << database_x_refs <<
            "Comments: "   << comments        <<
            "Kegg Terms: " << kegg_terms;
        return out.str();
    }

    bool is_empty() {
        return this->database_x_refs.empty() && this->comments.empty();
    }
};

class EntapDatabase {

public:

    // database typedefs
    typedef std::unordered_map<std::string, TaxEntry> tax_serial_map_t;
    typedef std::unordered_map<std::string, GoEntry> go_serial_map_t;
    typedef std::unordered_map<std::string, UniprotEntry> uniprot_serial_map_t;

    typedef enum {

        ENTAP_SERIALIZED=0, // Serialized database
        ENTAP_SQL,          // SQL database (uniprot mapping, tax data)

        ENTAP_MAX_TYPES

    } DATABASE_TYPE;

    typedef enum {

        ERR_DATA_OK=0,
        ERR_DATA_WRITE,
        ERR_DATA_READ,
        ERR_DATA_SQL_DUPLICATE,
        ERR_DATA_SQL_CREATE_DATABASE,
        ERR_DATA_SQL_TAX_CREATE_TABLE,      // 5
        ERR_DATA_SQL_CREATE_ENTRY,
        ERR_DATA_SQL_OPEN,
        ERR_DATA_FILE_EXISTS,
        ERR_DATA_TAX_CREATED,
        ERR_DATA_TAX_DOWNLOAD,              // 10
        ERR_DATA_FILE_DECOMPRESS,
        ERR_DATA_GO_DOWNLOAD,
        ERR_DATA_GO_DECOMPRESS,
        ERR_DATA_GO_ENTRY,
        ERR_DATA_SERIAL_FTP,                // 15
        ERR_DATA_SERIAL_DECOMPRESS,
        ERR_DATA_SQL_FTP,
        ERR_DATA_SQL_DECOMPRESS,
        ERR_DATA_FILE_MOVE,
        ERR_DATA_SERIALIZE_SAVE,            // 20
        ERR_DATA_SERIALIZE_READ,
        ERR_DATA_SERIAL_DUPLICATE,
        ERR_DATA_UNIPROT_DOWNLOAD,
        ERR_DATA_UNIPROT_DECOMPRESS,
        ERR_DATA_UNIPROT_ENTRY,             // 25
        ERR_DATA_UNIPROT_PARSE,
        ERR_DATA_UNIPROT_FILE,
        ERR_DATA_SQL_GO_CREATE_TABLE,
        ERR_DATA_SQL_UNIPROT_CREATE_TABLE,
        ERR_DATA_GO_PARSE,                  // 30
        ERR_DATA_TAXONOMY_PARSE,
        ERR_DATA_SET,
        ERR_DATA_INCOMPATIBLE_VER,
        ERR_DATA_GET_VERSION,
        ERR_DATA_DELETE_TABLE,
        ERR_DATA_CREATE_TABLE,

        ERR_DATA_MEM_ALLOC,
        ERR_DATA_UNHANDLED_TYPE,
        ERR_DATA_MAX=100

    } DATABASE_ERR;

    typedef enum {

        BOOST_TEXT_ARCHIVE=0,
        BOOST_BIN_ARCHIVE,
        CEREAL_BIN_ARCHIVE

    } SERIALIZATION_TYPE;

    struct EntapDatabaseStruct {
        tax_serial_map_t taxonomic_data;
        go_serial_map_t  gene_ontology_data;    // Accession - "GO:453232143"
        uniprot_serial_map_t uniprot_data;
        uint8 MAJOR_VERSION;
        uint8 MINOR_VERSION;

        EntapDatabaseStruct () {
            MAJOR_VERSION = 0;
            MINOR_VERSION = 0;
        };

#ifdef USE_BOOST
        friend class boost::serialization::access;
        template<typename Archive>
        void serialize(Archive & ar, const uint32 v) {
            ar&taxonomic_data;
            ar&gene_ontology_data;
            ar&uniprot_data;
            ar&MAJOR_VERSION;
            ar&MINOR_VERSION;
        }
#else
        // Use CEREAL for serialization
        template<class Archive>
        void serialize(Archive & archive)
        {
            archive(
                    taxonomic_data, gene_ontology_data, uniprot_data,
                    MAJOR_VERSION, MINOR_VERSION);
        }

#endif
    };

    // Node for NCBI taxonomy
    struct TaxonomyNode{
        std::string parent_id;
        std::string ncbi_id;
        std::string sci_name;
        vect_str_t  names;

        TaxonomyNode(std::string id);
    };

    EntapDatabase(FileSystem* fileSystem);
    ~EntapDatabase();
    bool set_database(DATABASE_TYPE type, std::string database_path);
    DATABASE_ERR download_database(DATABASE_TYPE, std::string&);
    DATABASE_ERR generate_database(DATABASE_TYPE, std::string&);
    std::string print_error_log();
    go_format_t format_go_delim(std::string terms, char delim);

    // Database accession routines
    TaxEntry get_tax_entry(std::string& species);
    GoEntry get_go_entry(std::string& go_id);
    UniprotEntry get_uniprot_entry(std::string& accession);

    bool is_uniprot_entry(std::string &sseqid, UniprotEntry &entry);

    // Database versioning
    bool is_valid_version();
    std::string get_current_version_str();
    std::string get_required_version_str();

#ifndef UNIT_TESTS
private:
#endif

    typedef enum {
        ENTAP_TAXONOMY,     // NCBI tax database
        ENTAP_GENE_ONTOLOGY,// GO database
        ENTAP_UNIPROT,      // UniProt mapping database
        ENTAP_VERSION,      // Version table used in SQL database only

        ENTAP_MAX_TABLES

    } DATABASE_TABLES;

    // Generation/download database routines
    DATABASE_ERR download_entap_sql(std::string&);
    DATABASE_ERR download_entap_serial(std::string&);
    DATABASE_ERR generate_entap_database(DATABASE_TYPE type, std::string& path);
    DATABASE_ERR generate_entap_tax(DATABASE_TYPE);
    DATABASE_ERR generate_entap_go(DATABASE_TYPE);
    DATABASE_ERR generate_entap_uniprot(DATABASE_TYPE);
    std::string  entap_tax_get_lineage(TaxonomyNode &,
                                       std::unordered_map<std::string, TaxonomyNode>&);

    bool create_sql_table(DATABASE_TABLES database_table);
    EntapDatabase::DATABASE_ERR create_database_type(DATABASE_TYPE type, std::string &path);
    EntapDatabase::DATABASE_ERR delete_database_table(DATABASE_TYPE type, DATABASE_TABLES table);


    bool sql_add_tax_entry(TaxEntry&);
    bool sql_add_go_entry(GoEntry&);
    bool add_uniprot_entry(DATABASE_TYPE type, UniprotEntry &entry);
    void set_err_msg(std::string msg, DATABASE_ERR code);
    bool set_database_versions(DATABASE_TYPE type);
    std::string get_uniprot_accession(std::string& sseqid);

    DATABASE_ERR serialize_database_save(SERIALIZATION_TYPE, std::string&);
    DATABASE_ERR serialize_database_read(SERIALIZATION_TYPE, std::string&);

    // FTP Paths
    const std::string FTP_GO_DATABASE =
            "http://archive.geneontology.org/latest-full/go_monthly-termdb-tables.tar.gz";
//    const std::string FTP_GO_DATABASE =
//            "http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz";
    const std::string FTP_NCBI_TAX_DUMP_TARGZ =
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
    const std::string FTP_ENTAP_DATABASE_SQL    =
            "https://treegenesdb.org/FTP/EnTAP/v0.10.8/databases/entap_database.db.gz";
    const std::string FTP_ENTAP_DATABASE_SERIAL =
            "https://treegenesdb.org/FTP/EnTAP/v0.10.8/databases/entap_database.bin.gz";
    const std::string FTP_UNIPROT_FLAT_FILE     =
            "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";

    const std::string ENTAP_DATABASE_SERIAL_GZ = "entap_database.bin.gz";
    const std::string ENTAP_DATABASE_SQL_GZ            = "entap_database.db.gz";

    // NCBI Taxonomy filenames
    const std::string NCBI_TAX_ROOT          = "root[Subtree]"; // Unused
    const std::string NCBI_TAX_DATABASE      = "taxonomy";      // Unused
    const std::string NCBI_TAX_DUMP_FILENAME = "taxdump.tar.gz";
    const std::string NCBI_TAX_DUMP_FTP_NAMES= "names.dmp";
    const std::string NCBI_TAX_DUMP_FTP_NODES= "nodes.dmp";
    const char        NCBI_TAX_DUMP_DELIM    = '\t';

    // NCBI Taxonomy dump columns
    const uint8 NCBI_TAX_DUMP_COL_NAME_CLASS   = 6; // Name type (scientific, authority...)
    const uint8 NCBI_TAX_DUMP_COL_ID           = 0;
    const uint8 NCBI_TAX_DUMP_COL_NAME         = 2; // Actual name
    const uint8 NCBI_TAX_DUMP_COL_PARENT       = 2; // Parent node in database
    const std::string NCBI_TAX_DUMP_SCIENTIFIC = "scientific name";

    // SQL Database Namings / Column numbers
    const std::string SQL_TABLE_NCBI_TAX_TITLE   = "TAXONOMY";
    const std::string SQL_COL_NCBI_TAX_TAXID     = "TAXID";
    const std::string SQL_COL_NCBI_TAX_LINEAGE   = "LINEAGE";
    const std::string SQL_COL_NCBI_TAX_NAME      = "TAXNAME"; // ex: homo sapiens

    const std::string SQL_TABLE_GO_TITLE         = "GENEONTOLOGY";
    const std::string SQL_TABLE_GO_COL_ID        = "GOID";
    const std::string SQL_TABLE_GO_COL_DESC      = "DESCRIPTION";
    const std::string SQL_TABLE_GO_COL_CATEGORY  = "CATEGORY";
    const std::string SQL_TABLE_GO_COL_LEVEL     = "LEVEL_INT"; // INTEGER
    const std::string SQL_TABLE_GO_COL_LEVEL_STR = "LEVEL_STR";

    const std::string SQL_TABLE_UNIPROT_TITLE    = "UNIPROT";
    const std::string SQL_TABLE_UNIPROT_COL_COMM = "COMMENTS";
    const std::string SQL_TABLE_UNIPROT_COL_XREF = "DATAXREFS";
    const std::string SQL_TABLE_UNIPROT_COL_ID   = "UNIPROTID";
    const std::string SQL_TABLE_UNIPROT_COL_KEGG = "KEGG";
    const std::string SQL_TABLE_UNIPROT_COL_GO   = "GOTERMS";

    const std::string SQL_TABLE_VERSION_TITLE    = "VERSION";
    const std::string SQL_TABLE_VERSION_COL_VER  = "VERSION";

    // Gene Ontology constants
    const std::string GO_BIOLOGICAL_LVL = "6679";
    const std::string GO_MOLECULAR_LVL  = "2892";
    const std::string GO_CELLULAR_LVL   = "311";
    const std::string GO_TERM_FILE      = "term.txt";
    const std::string GO_GRAPH_FILE     = "graph_path.txt";
//    const std::string GO_TERMDB_FILE    = "go_daily-termdb-tables.tar.gz";
    const std::string GO_TERMDB_FILE    = "go_monthly-termdb-tables.tar.gz";
//    const std::string GO_TERMDB_DIR     = "go_daily-termdb-tables/";
    const std::string GO_TERMDB_DIR     = "go_monthly-termdb-tables/";

    // UniProt mapping constants
    const std::string UNIPROT_DAT_FILE_GZ            = "uniprot_sprot.dat.gz";
    const std::string UNIPROT_DAT_FILE               = "uniprot_sprot.dat";

    // EnTAP database consts
    const SERIALIZATION_TYPE SERIALIZE_DEFAULT    = CEREAL_BIN_ARCHIVE;
    const uint8              SERIALIZE_MAJOR      = 2;
    const uint8              SERIALIZE_MINOR      = 0;
    const uint8              SQL_MAJOR            = 2;
    const uint8              SQL_MINOR            = 0;

    const uint8 STATUS_UPDATES = 5;     // Percentage of updates when downloading/configuring

    EntapDatabaseStruct *mpSerializedDatabase;
    FileSystem          *mpFileSystem;
    SQLDatabaseHelper   *mpDatabaseHelper;
    std::string          mTempDirectory;
    go_serial_map_t      mSqlGoHelper;    // Using to increase speeds for now, change later
    bool                 mUseSerial;
    std::string          mErrMsg;
    DATABASE_ERR         mErrCode;

    const std::string ENTAP_DATABASE_TYPES_STR[ENTAP_MAX_TYPES] {
            "EnTAP Serialized Database",
            "EnTAP SQL Database"
    };

    const std::string ENTAP_DATABASE_TABLES[ENTAP_MAX_TABLES] {
            "EnTAP NCBI Taxonomy Database",
            "EnTAP Gene Ontology Database",
            "EnTAP UniProt Swiss-Prot Database"
    };
};


#endif //ENTAP_ENTAPDATABASE_H
