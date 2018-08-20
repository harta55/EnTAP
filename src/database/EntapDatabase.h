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

#ifndef ENTAP_ENTAPDATABASE_H
#define ENTAP_ENTAPDATABASE_H

#include "../EntapGlobals.h"
#include "../EntapConfig.h"
#include "SQLDatabaseHelper.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "../FileSystem.h"


struct  GoEntry {
    std::string go_id;
    std::string level;
    std::string category;
    std::string term;
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const uint32 v) {
        ar&go_id;
        ar&level;
        ar&category;
        ar&term;
    }

    GoEntry() {
        go_id = "";
        level = "";
        category = "";
        term = "";
    }

    bool is_empty() {
        return term.empty();
    }
};

struct TaxEntry {
    std::string tax_id;
    std::string lineage;
    std::string tax_name;
    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive & ar, const uint32 v) {
        ar&tax_id;
        ar&lineage;
        ar&tax_name;
    }
    bool is_empty() {
        return this->tax_id.empty() && this->lineage.empty();
    }
    TaxEntry() {
        tax_id = "";
        lineage = "";
        tax_name = "";
    }
};

struct UniprotEntry {
    std::string database_x_refs;        // DR   OrthoDB; VOG090000I8; -.
    std::string comments;               // CC   -!- FUNCTION: Transcription activation. {ECO:0000305}.
    std::string uniprot_id;             // ID   001R_FRG3G              Reviewed;         256 AA.

    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive & ar, const uint32 v) {
        ar&uniprot_id;
        ar&database_x_refs;
        ar&comments;
    }

    UniprotEntry() {
        database_x_refs = "";
        comments = "";
        uniprot_id = "";
    }

    std::string print(void) {
        std::stringstream out;

        out <<
            "UniProt ID: " << uniprot_id <<
            "XRefs: "      << database_x_refs <<
            "Comments: "   << comments;
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
        ENTAP_TAXONOMY,     // NCBI tax database
        ENTAP_GENE_ONTOLOGY,// GO database
        ENTAP_UNIPROT,  // UniProt mapping database

        ENTAP_MAX_TYPES=5

    } DATABASE_TYPE;

    typedef enum {

        ERR_DATA_OK=0,
        ERR_DATA_WRITE,
        ERR_DATA_READ,
        ERR_DATA_SQL_DUPLICATE,
        ERR_DATA_SQL_CREATE_DATABASE,
        ERR_DATA_SQL_CREATE_TABLE,
        ERR_DATA_SQL_CREATE_ENTRY,
        ERR_DATA_SQL_OPEN,
        ERR_DATA_FILE_EXISTS,
        ERR_DATA_TAX_CREATED,
        ERR_DATA_TAX_DOWNLOAD,
        ERR_DATA_FILE_DECOMPRESS,
        ERR_DATA_GO_DOWNLOAD,
        ERR_DATA_GO_DECOMPRESS,
        ERR_DATA_GO_ENTRY,
        ERR_DATA_SERIAL_FTP,
        ERR_DATA_SERIAL_DECOMPRESS,
        ERR_DATA_SQL_FTP,
        ERR_DATA_SQL_DECOMPRESS,
        ERR_DATA_FILE_MOVE,
        ERR_DATA_SERIALIZE_SAVE,
        ERR_DATA_SERIALIZE_READ,
        ERR_DATA_SERIAL_DUPLICATE,
        ERR_DATA_UNIPROT_DOWNLOAD,
        ERR_DATA_UNIPROT_DECOMPRESS,
        ERR_DATA_UNIPROT_ENTRY,
        ERR_DATA_UNIPROT_PARSE,
        ERR_DATA_UNIPROT_FILE,

        ERR_DATA_MEM_ALLOC,
        ERR_DATA_MAX=100

    } DATABASE_ERR;

    typedef enum {

        BOOST_TEXT_ARCHIVE=0,
        BOOST_BIN_ARCHIVE

    } SERIALIZATION_TYPE;

    struct EntapDatabaseStruct {
        tax_serial_map_t taxonomic_data;
        go_serial_map_t  gene_ontology_data;
        uniprot_serial_map_t uniprot_data;

        friend class boost::serialization::access;
        template<typename Archive>
        void serialize(Archive & ar, const uint32 v) {
            ar&taxonomic_data;
            ar&gene_ontology_data;
            ar&uniprot_data;
        }
    };

    // Node for NCBI taxonomy
    struct TaxonomyNode{
        std::string parent_id;
        std::string ncbi_id;
        std::string sci_name;
        vect_str_t  names;

        TaxonomyNode(std::string id);
    };

    EntapDatabase(FileSystem*);
    ~EntapDatabase();
    bool set_database(DATABASE_TYPE, std::string);
    DATABASE_ERR download_database(DATABASE_TYPE, std::string&);
    DATABASE_ERR generate_database(DATABASE_TYPE, std::string&);
    std::string print_error_log(DATABASE_ERR err_code);

    // Database accession routines
    TaxEntry get_tax_entry(std::string& species);
    GoEntry get_go_entry(std::string& go_id);
    UniprotEntry get_uniprot_entry(std::string& accession);

    // Database accession routine (just making template)
//    template<class T>
//    T get_database_entry(std::string &accession) {
//        return T();
//    }


private:

    // Generation/download database routines
    DATABASE_ERR download_entap_sql(std::string&);
    DATABASE_ERR download_entap_serial(std::string&);
    DATABASE_ERR generate_entap_database(DATABASE_TYPE type, std::string& path);
    DATABASE_ERR generate_entap_tax(DATABASE_TYPE, std::string);
    DATABASE_ERR generate_entap_go(DATABASE_TYPE, std::string);
    DATABASE_ERR generate_entap_uniprot(DATABASE_TYPE, std::string);
    std::string  entap_tax_get_lineage(TaxonomyNode &,
                                       std::unordered_map<std::string, TaxonomyNode>&);
    bool sql_add_tax_entry(TaxEntry&);
    bool sql_add_go_entry(GoEntry&);
    bool create_sql_table(DATABASE_TYPE);
    bool add_uniprot_entry(DATABASE_TYPE type, UniprotEntry &entry);

    DATABASE_ERR serialize_database_save(SERIALIZATION_TYPE, std::string&);
    DATABASE_ERR serialize_database_read(SERIALIZATION_TYPE, std::string&);

    // FTP Paths
    const std::string FTP_GO_DATABASE =
            "http://archive.geneontology.org/latest-full/go_monthly-termdb-tables.tar.gz";
    const std::string FTP_NCBI_TAX_DUMP_TARGZ =
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
    const std::string FTP_ENTAP_DATABASE_SQL    = "https://treegenesdb.org/FTP/EnTAP/latest/databases/entap_database.db.gz";
    const std::string FTP_ENTAP_DATABASE_SERIAL = "https://treegenesdb.org/FTP/EnTAP/latest/databases/entap_database.bin.gz";
    const std::string FTP_UNIPROT_FLAT_FILE     = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";

    // NCBI Taxonomy filenames
    const std::string NCBI_TAX_ROOT          = "root[Subtree]"; // Unused
    const std::string NCBI_TAX_DATABASE      = "taxonomy";      // Unused
    const std::string NCBI_TAX_DUMP_FILENAME = "taxdump.tar.gz";
    const std::string NCBI_TAX_DUMP_FTP_NAMES= "names.dmp";
    const std::string NCBI_TAX_DUMP_FTP_NODES= "nodes.dmp";
    const char        NCBI_TAX_DUMP_DELIM    = '\t';

    // NCBI Taxonomy dump columns
    const uint8 NCBI_TAX_DUMP_COL_NAME_CLASS = 6; // Name type (scientific, authority...)
    const uint8 NCBI_TAX_DUMP_COL_ID         = 0;
    const uint8 NCBI_TAX_DUMP_COL_NAME       = 2; // Actual name
    const uint8 NCBI_TAX_DUMP_COL_PARENT     = 2; // Parent node in database
    const std::string NCBI_TAX_DUMP_SCIENTIFIC = "scientific name";

    // SQL Database Namings / Column numbers
    const std::string SQL_TABLE_NCBI_TAX_TITLE = "TAXONOMY";
    const std::string SQL_COL_NCBI_TAX_TAXID   = "TAXID";
    const std::string SQL_COL_NCBI_TAX_LINEAGE = "LINEAGE";
    const std::string SQL_COL_NCBI_TAX_NAME    = "TAXNAME"; // ex: homo sapiens
    const std::string SQL_TABLE_GO_TITLE       = "GENEONTOLOGY";
    const std::string SQL_TABLE_GO_COL_ID      = "GOID";
    const std::string SQL_TABLE_GO_COL_DESC    = "DESCRIPTION";
    const std::string SQL_TABLE_GO_COL_CATEGORY= "CATEGORY";
    const std::string SQL_TABLE_GO_COL_LEVEL   = "LEVEL";
    const std::string SQL_TABLE_UNIPROT_TITLE  = "UNIPROT";
    const std::string SQL_TABLE_UNIPROT_COL_COMM= "COMMENTS";
    const std::string SQL_TABLE_UNIPROT_COL_XREF= "DATAXREFS";
    const std::string SQL_TABLE_UNIPROT_COL_ID  = "UNIPROTID";

    // Gene Ontology constants
    const std::string GO_BIOLOGICAL_LVL = "6679";
    const std::string GO_MOLECULAR_LVL  = "2892";
    const std::string GO_CELLULAR_LVL   = "311";
    const std::string GO_TERM_FILE      = "term.txt";
    const std::string GO_GRAPH_FILE     = "graph_path.txt";
    const std::string GO_TERMDB_FILE    = "go_monthly-termdb-tables.tar.gz";
    const std::string GO_TERMDB_DIR     = "go_monthly-termdb-tables/";

    // UniProt mapping constants
    const std::string UNIPROT_DAT_FILE_GZ            = "uniprot_sprot.dat.gz";
    const std::string UNIPROT_DAT_FILE               = "uniprot_sprot.dat";
    const uint8       UNIPROT_DAT_TAG_LEN            = 2;   // Length of tags
    const uint8       UNIPROT_DAT_TAG_DATA_POS       = 5;   // Position data starts
    // ID   001R_FRG3G              Reviewed;         256 AA.
    const std::string UNIPROT_DAT_TAG_ID             = "ID";
    // DR   SwissPalm; Q6GZX4; -.
    const std::string UNIPROT_DAT_TAG_DATABASE_X_REF = "DR";
    // CC   -!- FUNCTION: Transcription activation. {ECO:0000305}.
    const std::string UNIPROT_DAT_TAG_COMMENT        = "CC";
    // Entries are split by this (this is on the last line of file)
    const std::string UNIPROT_DAT_TAG_NEXT_ENTRY     = "//";

    // EnTAP database consts
    const SERIALIZATION_TYPE SERIALIZE_DEFAULT    = BOOST_TEXT_ARCHIVE;

    const uint8 STATUS_UPDATES = 5;     // Percentage of updates when downloading/configuring

    EntapDatabaseStruct *_pSerializedDatabase;
    FileSystem          *_pFilesystem;
    SQLDatabaseHelper   *_pDatabaseHelper;
    std::string          _temp_directory;
    go_serial_map_t      _sql_go_helper;    // Using to increase speeds for now, change later
    bool                 _use_serial;
    std::string          _err_msg;

    const std::string ENTAP_DATABASE_TYPES_STR[ENTAP_MAX_TYPES] {
            "EnTAP Serialized Database",
            "EnTAP SQL Database",
            "EnTAP NCBI Taxonomy Database",
            "EnTAP Gene Ontology Database",
            "EnTAP UniProt Swiss-Prot Database"
    };
};


#endif //ENTAP_ENTAPDATABASE_H
