//
// Created by harta on 3/7/18.
//

#ifndef ENTAP_ENTAPDATABASE_H
#define ENTAP_ENTAPDATABASE_H

#include "../common.h"
#include "../EntapConfig.h"

class EntapDatabase {

public:
    struct EntapDatabaseStruct {
        tax_serial_map_t taxonomic_data;
        go_serial_map_t  gene_ontology_data;
    };

    EntapDatabase(std::string&);

private:
    EntapDatabaseStruct *DATABASE;
};


#endif //ENTAP_ENTAPDATABASE_H
