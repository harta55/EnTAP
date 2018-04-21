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
        DATA_SQL_CREATE,
        DATA_FILE_EXISTS,
        DATA_TAX_CREATION

    } DATABASE_ERR;

    struct EntapDatabaseStruct {
        tax_serial_map_t taxonomic_data;
        go_serial_map_t  gene_ontology_data;
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


    const std::string FTP_ENTAP_DATABASE_SQL = "temp/ftp.com";
    const std::string NCBI_TAX_ROOT          = "root[Subtree]"; // Unused
    const std::string NCBI_TAX_DATABASE      = "taxonomy";      // Unused
    const std::string NCBI_TAX_DUMP_FTP_TARGZ= "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
    const std::string NCBI_TAX_DUMP_FTP_NAMES= "names.dmp";
    const std::string NCBI_TAX_DUMP_FTP_NODES= "nodes.dmp";

    EntapDatabaseStruct *_pDATABASE;
    FileSystem          *_pFilesystem;
    SQLDatabaseHelper   *_pDatabaseHelper;
    std::string          _temp_directory;


};


#endif //ENTAP_ENTAPDATABASE_H
