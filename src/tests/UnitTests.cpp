/*
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
#include "../config.h"
#ifdef UNIT_TESTS
#include "UnitTests.h"
#include "../database/EntapDatabase.h"
#include "../QueryData.h"
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

    // Test suite of QueryData tests
//    TestQueryData();

}

void UnitTests::TestEntapDatabase() {

    auto *entapDatabase = new EntapDatabase(mpFileSystem);

    try {
        UTEntapDatabase_00(entapDatabase, EntapDatabase::ENTAP_SQL);      // SQL Generate
        UTEntapDatabase_01(entapDatabase, EntapDatabase::ENTAP_SQL);      // SQL Gene Ontology
        UTEntapDatabase_02(entapDatabase, EntapDatabase::ENTAP_SQL);      // SQL Taxonomy
        UTEntapDatabase_03(entapDatabase, EntapDatabase::ENTAP_SQL);      // SQL UniProt

        UTEntapDatabase_00(entapDatabase, EntapDatabase::ENTAP_SERIALIZED);      // Serialized Generate
        UTEntapDatabase_01(entapDatabase, EntapDatabase::ENTAP_SERIALIZED);      // Serialized Gene Ontology
        UTEntapDatabase_02(entapDatabase, EntapDatabase::ENTAP_SERIALIZED);      // Serialized Taxonomy
        UTEntapDatabase_03(entapDatabase, EntapDatabase::ENTAP_SERIALIZED);      // Serialized UniProt

    } catch (...) {
        delete entapDatabase;
    }
}

void UnitTests::UTEntapDatabase_00(EntapDatabase *entapDatabase, EntapDatabase::DATABASE_TYPE type) {
    EntapDatabase::DATABASE_ERR databaseErr;
    std::string database_path = PATHS(mpFileSystem->get_temp_outdir(), "database_temp.db");

    if (mpFileSystem->file_exists(database_path)) {
        mpFileSystem->delete_file(database_path);
    }

    databaseErr = entapDatabase->create_database_type(type, database_path);
    assert(databaseErr == EntapDatabase::ERR_DATA_OK);
    std::cout << "UT_EntapDatabase_00 COMPLETE" << std::endl;
}

// SQL Gene Ontology Table tests
void UnitTests::UTEntapDatabase_01(EntapDatabase *entapDatabase, EntapDatabase::DATABASE_TYPE type) {
    EntapDatabase::DATABASE_ERR databaseErr;
    GoEntry goEntry;

    databaseErr = entapDatabase->generate_entap_go(type);
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

    FS_dprint("UT_EntapDatabase_01 COMPLETE - Type: " + std::to_string(type));
    std::cout << "UT_EntapDatabase_01 COMPLETE" << std::endl;
}

// SQL Taxonomy Table tests
void UnitTests::UTEntapDatabase_02(EntapDatabase *entapDatabase, EntapDatabase::DATABASE_TYPE type) {
    EntapDatabase::DATABASE_ERR databaseErr;
    TaxEntry taxEntry;
    std::string tax_test;

    databaseErr = entapDatabase->generate_entap_tax(type);
    assert(databaseErr == EntapDatabase::ERR_DATA_OK);

    tax_test = "pinus";
    taxEntry = entapDatabase->get_tax_entry(tax_test);
    assert(taxEntry.tax_name == "pinus");
    assert(taxEntry.lineage == "pinus;pinaceae;pinales;pinidae;pinopsida;acrogymnospermae;"
                               "spermatophyta;euphyllophyta;tracheophyta;embryophyta;streptophytina;"
                               "streptophyta;viridiplantae;eukaryota;cellular organisms;root");
    assert(taxEntry.tax_id == "3337");

    tax_test = "homo sapiens";
    taxEntry = entapDatabase->get_tax_entry(tax_test);
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
    std::cout << "UT_EntapDatabase_02 - Type: "+ std::to_string(type) << std::endl;
}

// SQL UniProt
void UnitTests::UTEntapDatabase_03(EntapDatabase *entapDatabase, EntapDatabase::DATABASE_TYPE type) {
    EntapDatabase::DATABASE_ERR databaseErr;
    UniprotEntry uniprotEntry;
    std::string uniprot_test;

    databaseErr = entapDatabase->generate_entap_uniprot(type);
    assert(databaseErr == EntapDatabase::ERR_DATA_OK);

    uniprot_test = "MK67I_BOVIN";
    uniprotEntry = entapDatabase->get_uniprot_entry(uniprot_test);
    assert(uniprotEntry.uniprot_id == "MK67I_BOVIN");
    assert(uniprotEntry.go_terms.size() == 6);
    for (const GoEntry &entry : uniprotEntry.go_terms) {
        assert(entry.go_id == "GO:0003723" || entry.go_id == "GO:0000463" || entry.go_id == "GO:0005654" ||
               entry.go_id == "GO:0005730" || entry.go_id == "GO:0005737" || entry.go_id == "GO:0000794");
    }

    uniprot_test = "ANFY1_MOUSE";
    uniprotEntry = entapDatabase->get_uniprot_entry(uniprot_test);
    assert(uniprotEntry.uniprot_id == "ANFY1_MOUSE");
    assert(uniprotEntry.go_terms.size() == 17);
    for (const GoEntry &entry : uniprotEntry.go_terms) {
        assert(entry.go_id == "GO:0005829" || entry.go_id == "GO:0005769" || entry.go_id == "GO:0005768" ||
               entry.go_id == "GO:0010008" || entry.go_id == "GO:0043231" || entry.go_id == "GO:0044354" ||
               entry.go_id == "GO:0016020" || entry.go_id == "GO:0030904" || entry.go_id == "GO:0046872" ||
               entry.go_id == "GO:1901981" || entry.go_id == "GO:0017137" || entry.go_id == "GO:0006897" ||
               entry.go_id == "GO:0016197" || entry.go_id == "GO:0034058" || entry.go_id == "GO:0090160" ||
               entry.go_id == "GO:0048549" || entry.go_id == "GO:0042147");
    }


//    databaseErr = entapDatabase->delete_database_table(EntapDatabase::ENTAP_SQL, EntapDatabase::ENTAP_GENE_ONTOLOGY);
//    assert(databaseErr == EntapDatabase::ERR_DATA_OK);

    FS_dprint("UT_EntapDatabase_01 COMPLETE");
    std::cout << "UT_EntapDatabase_01 COMPLETE:  - Type: " + std::to_string(type) << std::endl;
}

void UnitTests::TestQueryData() {
    QueryData queryData = QueryData();
    std::string base_path = PATHS(mRootDirectory, "test_alignment");
    std::vector<ENTAP_HEADERS> header = {
            ENTAP_HEADER_ONT_EGG_SEED_ORTHO,
            ENTAP_HEADER_ONT_EGG_SEED_EVAL,
            ENTAP_HEADER_ONT_EGG_SEED_SCORE,
            ENTAP_HEADER_ONT_EGG_PRED_GENE,
            ENTAP_HEADER_ONT_EGG_TAX_SCOPE_READABLE,
            ENTAP_HEADER_ONT_EGG_TAX_SCOPE_MAX,
            ENTAP_HEADER_ONT_EGG_MEMBER_OGS,
            ENTAP_HEADER_ONT_EGG_DESC,
            ENTAP_HEADER_ONT_EGG_KEGG,
            ENTAP_HEADER_ONT_EGG_BIGG,
            ENTAP_HEADER_ONT_EGG_GO_BIO,
            ENTAP_HEADER_ONT_EGG_GO_CELL,
            ENTAP_HEADER_ONT_EGG_GO_MOLE,
            ENTAP_HEADER_ONT_EGG_PROTEIN
    };

    vect_uint16_t  go_levels = {0,1};
    std::vector<FileSystem::ENT_FILE_TYPES> types = std::vector<FileSystem::ENT_FILE_TYPES>{FileSystem::ENT_FILE_DELIM_TSV,
                                                            FileSystem::ENT_FILE_FASTA_FAA,
                                                            FileSystem::ENT_FILE_FASTA_FNN,
                                                            FileSystem::ENT_FILE_GENE_ENRICH_EFF_LEN,
                                                            FileSystem::ENT_FILE_GENE_ENRICH_GO_TERM};
    queryData.start_alignment_files(base_path, header,go_levels,types);
    queryData.end_alignment_files(base_path);
}

void UnitTests::TestEggnogDatabase() {

}

#endif