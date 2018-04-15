//
// Created by harta on 3/7/18.
//

#ifndef ENTAP_ENTAPDATABASE_H
#define ENTAP_ENTAPDATABASE_H

#include "../common.h"
#include "../EntapConfig.h"

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
        DATA_WRITE

    } DATABASE_ERR;

    struct EntapDatabaseStruct {
        tax_serial_map_t taxonomic_data;
        go_serial_map_t  gene_ontology_data;
    };

    EntapDatabase();
    void set_database(DATABASE_TYPE);

private:
    EntapDatabaseStruct *DATABASE;
};


#endif //ENTAP_ENTAPDATABASE_H
