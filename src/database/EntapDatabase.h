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

#include "../common.h"
#include "../EntapConfig.h"
#include "SQLDatabaseHelper.h"

class EntapDatabase {

public:

    typedef enum {

        ENTAP_SQL=0,        // SQL database (uniprot mapping, tax data)
        ENTAP_TAXONOMY,     // NCBI tax database
        ENTAP_GENE_ONTOLOGY,// GO database
        ENTAP_UNIPROT_MAP,  // UniProt mapping database
        ENTAP_SERIALIZED,   // Serialized database (uniprot, tax, go)
        EGGNOG_SQL,         // EggNOG SQL database
        EGGNOG_DIAMOND,     // EggNOG DIAMOND database
        EGGNOG_FASTA        // EggNOG fasta (non-diamond)

    } DATABASE_TYPE;

    typedef enum {

        ERR_DATA_OK=0,
        ERR_DATA_WRITE,
        ERR_DATA_READ,
        ERR_DATA_SQL_DUPLICATE,
        ERR_DATA_SQL_CREATE_DATABASE,
        ERR_DATA_SQL_CREATE_TABLE,
        ERR_DATA_SQL_CREATE_ENTRY,
        ERR_DATA_FILE_EXISTS,
        ERR_DATA_TAX_CREATED,
        ERR_DATA_TAX_DOWNLOAD,
        ERR_DATA_FILE_DECOMPRESS,
        ERR_DATA_GO_DOWNLOAD,
        ERR_DATA_GO_DECOMPRESS,
        ERR_DATA_GO_ENTRY,
        ERR_DATA_SERIAL_FTP,
        ERR_DATA_SERIAL_DECOMPRESS,
        ERR_DATA_FILE_MOVE

    } DATABASE_ERR;

    struct EntapDatabaseStruct {
        tax_serial_map_t taxonomic_data;
        go_serial_map_t  gene_ontology_data;
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
    DATABASE_ERR set_database(DATABASE_TYPE);
    DATABASE_ERR download_database(DATABASE_TYPE, std::string&);
    DATABASE_ERR generate_database(DATABASE_TYPE, std::string&);

    // Database accession routine (just making template)
    template<class T>
    T get_database_entry(DATABASE_TYPE type, std::string &accession) {
        bool use_serial;

        if (accession.empty()) return T();

        // Check for presence of serialized first
        if (_pSerializedDatabase != nullptr) {
            use_serial = true;
        } else if (_pDatabaseHelper != nullptr) {
            use_serial = false;
        } else {
            return T();
        }

        switch (type) {
            case ENTAP_TAXONOMY:
                return get_tax_entry(accession, use_serial);
            case ENTAP_GENE_ONTOLOGY:
                return get_go_entry(accession, use_serial);
            default:
                return NULL;
        }
    }

private:
    // Generation/download database routines
    DATABASE_ERR download_entap_sql(std::string&);
    DATABASE_ERR download_entap_serial(std::string&);
    DATABASE_ERR generate_entap_sql(std::string&);
    DATABASE_ERR generate_entap_serial(std::string&);
    DATABASE_ERR generate_entap_tax(DATABASE_TYPE, std::string);
    DATABASE_ERR generate_entap_go(DATABASE_TYPE, std::string);
    std::string  entap_tax_get_lineage(TaxonomyNode &,
                                       std::unordered_map<std::string, TaxonomyNode>&);
    bool sql_add_tax_entry(TaxEntry&);
    bool sql_add_go_entry(GoEntry&);
    bool create_sql_table(DATABASE_TYPE);

    // Database accession routines
    TaxEntry get_tax_entry(std::string& species, bool use_serial);
    GoEntry get_go_entry(std::string& go_id, bool use_serial);

    // FTP Paths
    const std::string FTP_GO_DATABASE =
            "http://archive.geneontology.org/latest-full/go_monthly-termdb-tables.tar.gz";
    const std::string FTP_NCBI_TAX_DUMP_TARGZ =
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
    const std::string FTP_ENTAP_DATABASE_SQL = "temp/ftp";
    const std::string FTP_ENTAP_DATABASE_SERIAL = "temp.tar.gz";

    // NCBI Taxonomy filenames
    const std::string NCBI_TAX_ROOT          = "root[Subtree]"; // Unused
    const std::string NCBI_TAX_DATABASE      = "taxonomy";      // Unused
    const std::string NCBI_TAX_DUMP_FILENAME = "taxdump.tar.gz";
    const std::string NCBI_TAX_DUMP_FTP_NAMES= "names.dmp";
    const std::string NCBI_TAX_DUMP_FTP_NODES= "nodes.dmp";
    const char        NCBI_TAX_DUMP_DELIM    = '\t';

    // NCBI Taxonomy dump columns
    const int NCBI_TAX_DUMP_COL_NAME_CLASS = 6; // Name type (scientific, authority...)
    const int NCBI_TAX_DUMP_COL_ID         = 0;
    const int NCBI_TAX_DUMP_COL_NAME       = 2; // Actual name
    const int NCBI_TAX_DUMP_COL_PARENT     = 2; // Parent node in database
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

    // Gene Ontology constants
    const std::string GO_BIOLOGICAL_LVL = "6679";
    const std::string GO_MOLECULAR_LVL  = "2892";
    const std::string GO_CELLULAR_LVL   = "311";
    const std::string GO_TERM_FILE      = "term.txt";
    const std::string GO_GRAPH_FILE     = "graph_path.txt";
    const std::string GO_TERMDB_FILE    = "go_monthly-termdb-tables.tar.gz";
    const std::string GO_TERMDB_DIR     = "go_monthly-termdb-tables/";

    // EnTAP database consts
    const std::string ENTAP_SERIAL_FILENAME_TARGZ = "entap_database.bin.tar.gz"; // temp name
    const std::string ENTAP_SERIAL_FILENAME       = "entap_database.bin";

    const uint8 STATUS_UPDATES = 5;     // Percentage of updates when downloading/configuring

    EntapDatabaseStruct *_pSerializedDatabase;
    FileSystem          *_pFilesystem;
    SQLDatabaseHelper   *_pDatabaseHelper;
    std::string          _temp_directory;


};


#endif //ENTAP_ENTAPDATABASE_H
