/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2020, Alexander Hart, Dr. Jill Wegrzyn
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
#include "../config.h"
#ifdef UNIT_TESTS
#include "UnitTests.h"
#include "../database/EntapDatabase.h"
#include <assert.h>

UnitTests::UnitTests() {
    mpFileSystem = new FileSystem();
    mRootDirectory = PATHS(mpFileSystem->get_root_path(), "Testing");
    mpFileSystem->set_root_dir(mRootDirectory);
}

UnitTests::~UnitTests() {
    delete mpFileSystem;
}

void UnitTests::execute_tests() {

    // Test suite of EnTAP database tests
    TestEntapDatabase();

}

void UnitTests::TestEntapDatabase() {

    auto *entapDatabase = new EntapDatabase(mpFileSystem);

    try {
        UTEntapDatabase_00(entapDatabase);      // SQL Generate
//        UTEntapDatabase_01(entapDatabase);      // SQL Gene Ontology
//        UTEntapDatabase_02(entapDatabase);      // SQL Taxonomy
        UTEntapDatabase_03(entapDatabase);      // SQL Taxonomy

    } catch (...) {
        delete entapDatabase;
    }
}

void UnitTests::UTEntapDatabase_00(EntapDatabase *entapDatabase) {
    EntapDatabase::DATABASE_ERR databaseErr;
    std::string database_path = PATHS(mpFileSystem->get_temp_outdir(), "database_temp.db");

    databaseErr = entapDatabase->create_database_type(EntapDatabase::ENTAP_SQL, database_path);
    assert(databaseErr == EntapDatabase::ERR_DATA_OK);
}

// SQL Gene Ontology Table tests
void UnitTests::UTEntapDatabase_01(EntapDatabase *entapDatabase) {
    EntapDatabase::DATABASE_ERR databaseErr;
    GoEntry goEntry;

    databaseErr = entapDatabase->generate_entap_go(EntapDatabase::ENTAP_SQL);
    assert(databaseErr == EntapDatabase::ERR_DATA_OK);

    std::string go_test;

    go_test = "GO:0000166";
    goEntry = entapDatabase->get_go_entry(go_test);
    assert(goEntry.category == "molecular_function");
    assert(goEntry.go_id == "GO:0000166");
    assert(goEntry.term == "nucleotide binding");

    go_test = "GO:0007275";
    goEntry = entapDatabase->get_go_entry(go_test);
    assert(goEntry.category == "biological_process");
    assert(goEntry.go_id == "GO:0007275");
    assert(goEntry.term == "multicellular organism development");

    go_test = "GO:0019717";
    goEntry = entapDatabase->get_go_entry(go_test);
    assert(goEntry.category == "cellular_component");
    assert(goEntry.go_id == "GO:0019717");
    assert(goEntry.term == "obsolete synaptosome");

//    databaseErr = entapDatabase->delete_database_table(EntapDatabase::ENTAP_SQL, EntapDatabase::ENTAP_GENE_ONTOLOGY);
//    assert(databaseErr == EntapDatabase::ERR_DATA_OK);

    FS_dprint("UT_EntapDatabase_01 COMPLETE");
}

// SQL Taxonomy Table tests
void UnitTests::UTEntapDatabase_02(EntapDatabase *entapDatabase) {
    EntapDatabase::DATABASE_ERR databaseErr;
    TaxEntry taxEntry;
    std::string tax_test;

    databaseErr = entapDatabase->generate_entap_tax(EntapDatabase::ENTAP_SQL);
    assert(databaseErr == EntapDatabase::ERR_DATA_OK);

    tax_test = "pinus";
    taxEntry = entapDatabase->get_tax_entry(tax_test);
    assert(taxEntry.tax_name == "pinus");
    assert(taxEntry.lineage == "pinus;pinaceae;pinales;pinidae;pinopsida;acrogymnospermae;"
                               "spermatophyta;euphyllophyta;tracheophyta;embryophyta;streptophytina;"
                               "streptophyta;viridiplantae;eukaryota;cellular organisms;root");
    assert(taxEntry.tax_id == "3337");

    tax_test = "homo sapiens";
    assert(taxEntry.tax_name == "homo sapiens");
    assert(taxEntry.lineage == "homo;homininae;hominidae;hominoidea;catarrhini;simiiformes;"
                               "haplorrhini;primates;euarchontoglires;boreoeutheria;eutheria;theria;"
                               "mammalia;amniota;tetrapoda;dipnotetrapodmorpha;sarcopterygii;"
                               "euteleostomi;teleostomi;gnathostomata;vertebrata;craniata;"
                               "deuterostomia;bilateria;eumetazoa;metazoa;opisthokonta;eukaryota;"
                               "cellular organisms;root");
    assert(taxEntry.tax_id == "9606");

//    databaseErr = entapDatabase->delete_database_table(EntapDatabase::ENTAP_SQL, EntapDatabase::ENTAP_TAXONOMY);
//    assert(databaseErr == EntapDatabase::ERR_DATA_OK);

    FS_dprint("UT_EntapDatabase_02 COMPLETE");
}

// SQL UniProt
void UnitTests::UTEntapDatabase_03(EntapDatabase *entapDatabase) {
    EntapDatabase::DATABASE_ERR databaseErr;
    UniprotEntry uniprotEntry;
    std::string uniprot_test;

    databaseErr = entapDatabase->generate_entap_uniprot(EntapDatabase::ENTAP_SQL);
    assert(databaseErr == EntapDatabase::ERR_DATA_OK);

    uniprot_test = "Q3SZM1";
    uniprotEntry = entapDatabase->get_uniprot_entry(uniprot_test);
    assert(uniprotEntry.uniprot_id == "Q3SZM1");
    assert(uniprotEntry.go_terms.begin()->go_id == "GO:0003723");
    assert(uniprotEntry.go_terms.begin()++->go_id == "GO:0000463");

//    databaseErr = entapDatabase->delete_database_table(EntapDatabase::ENTAP_SQL, EntapDatabase::ENTAP_GENE_ONTOLOGY);
//    assert(databaseErr == EntapDatabase::ERR_DATA_OK);

    FS_dprint("UT_EntapDatabase_01 COMPLETE");
}
#endif