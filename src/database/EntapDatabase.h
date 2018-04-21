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
        ENTAP_GENE_ONTOLOGY,// GO database (serial)
        EGGNOG_SQL,         // EggNOG SQL database
        EGGNOG_DIAMOND,     // EggNOG DIAMOND database
    } DATABASE_TYPE;

    typedef enum {

        DATA_OK=0,
        DATA_WRITE,
        DATA_READ,
        DATA_SQL_DUPLICATE,
        DATA_SQL_CREATE_DATABASE,
        DATA_SQL_CREATE_TABLE,
        DATA_SQL_CREATE_ENTRY,
        DATA_FILE_EXISTS,
        DATA_TAX_CREATION

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

private:
    DATABASE_ERR download_entap_sql(std::string&);
    DATABASE_ERR generate_entap_sql(std::string&);
    DATABASE_ERR generate_entap_tax(DATABASE_TYPE, std::string);
    DATABASE_ERR generate_entap_go(DATABASE_TYPE, std::string);
    std::string  entap_tax_get_lineage(TaxonomyNode &,
                                       std::unordered_map<std::string, TaxonomyNode>&);
    bool sql_add_tax_entry(TaxEntry&);
    bool create_sql_table(DATABASE_TYPE);

    const std::string FTP_ENTAP_DATABASE_SQL = "temp/ftp.com";
    const std::string NCBI_TAX_ROOT          = "root[Subtree]"; // Unused
    const std::string NCBI_TAX_DATABASE      = "taxonomy";      // Unused
    const std::string NCBI_TAX_DUMP_FTP_TARGZ= "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
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

    // SQL Database Namings
    const std::string SQL_TABLE_NCBI_TAX_TITLE = "TAXONOMY";
    const std::string SQL_COL_NCBI_TAX_TAXID   = "TAXID";
    const std::string SQL_COL_NCBI_TAX_LINEAGE = "LINEAGE";
    const std::string SQL_COL_NCBI_TAX_NAME    = "TAXNAME";
    const std::string SQL_TABLE_GO_TITLE       = "GENEONTOLOGY";

    EntapDatabaseStruct *_pSerializedDatabase;
    FileSystem          *_pFilesystem;
    SQLDatabaseHelper   *_pDatabaseHelper;
    std::string          _temp_directory;


};


#endif //ENTAP_ENTAPDATABASE_H
