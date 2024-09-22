/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2024, Alexander Hart, Dr. Jill Wegrzyn
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
#include "../database/NCBIDatabase.h"


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
    // TestEntapDatabase();

    // Test suite of QueryData tests
//    TestQueryData();

    // Test suite of NCBI Entrez tests
    // TestNCBIEntrez();

    // Test suite of NCBI Database tests
    TestNCBIDatabase();

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

        UTEntapDatabase_04(entapDatabase);      // Check NCBI accession routines

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

    FS_dprint("UT_EntapDatabase_03 COMPLETE");
    std::cout << "UT_EntapDatabase_03 COMPLETE:  - Type: " + std::to_string(type) << std::endl;
}

void UnitTests::UTEntapDatabase_04(EntapDatabase *entapDatabase) {
    std::string taxon;

    taxon = "homo sapiens";
    assert(entapDatabase->is_ncbi_tax_entry(taxon) == true);

    taxon = "homo_sapiens";
    assert(entapDatabase->is_ncbi_tax_entry(taxon) == true);

    taxon = "hfdasfa";
    assert(entapDatabase->is_ncbi_tax_entry(taxon) == false);

    FS_dprint("UT_EntapDatabase_04 COMPLETE");
    std::cout << "UT_EntapDatabase_04 COMPLETE" << std::endl;
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

void UnitTests::TestNCBIEntrez() {

    UTNCBIEntrez_00();

}

// Test Entrez fetch
void UnitTests::UTNCBIEntrez_00() {
    NCBIEntrez ncbi_entrez = NCBIEntrez(mpFileSystem);

    FS_dprint("UT_NCBIEntrez_00 START");
    std::cout << "UT_NCBIEntrez_00 START" << std::endl;

    // Test 1, pulling gene id from protein database for a single NCBI ID
    std::string test_ncbi_id = "XP_014245616";
    std::string test_expected_geneid = "106664428";
    NCBIEntrez::EntrezResults entrez_results;
    NCBIEntrez::EntrezInput entrez_input;
    entrez_input.data_types = {NCBIEntrez::ENTREZ_DATA_GENEID};
    entrez_input.uid_list = {test_ncbi_id};
    entrez_input.database = NCBIEntrez::NCBI_DATABASE_PROTEIN;
    entrez_input.rettype = NCBIEntrez::NCBI_ENTREZ_RETTYPE_GP;

    ncbi_entrez.entrez_fetch(entrez_input, entrez_results);
    assert(!entrez_results.entrez_results.empty());
    assert(entrez_results.entrez_results.find(test_ncbi_id) != entrez_results.entrez_results.end());
    assert(entrez_results.entrez_results.at(test_ncbi_id).geneid == test_expected_geneid);
    FS_dprint("\tUT_NCBIEntrez_00 - CASE_01 PASS");
    std::cout << "\tUT_NCBIEntrez_00 - CASE_01 PASS" << std::endl;


    // Test 2, pulling gene id from protein database for multiple NCBI ID
    //  including NCBI versions (XXXX.1)
    entrez_input = {};
    entrez_results = {};
    entrez_input.data_types = {NCBIEntrez::ENTREZ_DATA_GENEID};
    entrez_input.uid_list = {"XP_014245616", "XP_020482136.1", "XP_026476231.1"};
    entrez_input.database = NCBIEntrez::NCBI_DATABASE_PROTEIN;
    entrez_input.rettype = NCBIEntrez::NCBI_ENTREZ_RETTYPE_GP;

    ncbi_entrez.entrez_fetch(entrez_input, entrez_results);
    assert(!entrez_results.entrez_results.empty());
    assert(entrez_results.entrez_results.find("XP_014245616") != entrez_results.entrez_results.end());
    assert(entrez_results.entrez_results.at("XP_014245616").geneid == "106664428");
    assert(entrez_results.entrez_results.find("XP_020482136.1") != entrez_results.entrez_results.end());
    assert(entrez_results.entrez_results.at("XP_020482136.1").geneid == "109976342");
    assert(entrez_results.entrez_results.find("XP_026476231.1") != entrez_results.entrez_results.end());
    assert(entrez_results.entrez_results.at("XP_026476231.1").geneid == "113381705");
    FS_dprint("\tUT_NCBIEntrez_00 - CASE_02 PASS");
    std::cout << "\tUT_NCBIEntrez_00 - CASE_02 PASS" << std::endl;

    // Test 3, entrez_fetch returns FALSE if no data retrieved
    entrez_input = {};
    entrez_results = {};
    entrez_input.data_types = {NCBIEntrez::ENTREZ_DATA_GENEID};
    entrez_input.uid_list = {"XP_014245616safdsafs"};
    entrez_input.database = NCBIEntrez::NCBI_DATABASE_PROTEIN;
    entrez_input.rettype = NCBIEntrez::NCBI_ENTREZ_RETTYPE_GP;

    assert(false == ncbi_entrez.entrez_fetch(entrez_input, entrez_results));
    FS_dprint("\tUT_NCBIEntrez_00 - CASE_03 PASS");
    std::cout << "\tUT_NCBIEntrez_00 - CASE_03 PASS" << std::endl;

    // Test 4, entrez_fetch returns FALSE with no UID list
    entrez_input = {};
    entrez_results = {};
    entrez_input.data_types = {NCBIEntrez::ENTREZ_DATA_GENEID};
    entrez_input.uid_list = {};
    entrez_input.database = NCBIEntrez::NCBI_DATABASE_PROTEIN;
    entrez_input.rettype = NCBIEntrez::NCBI_ENTREZ_RETTYPE_GP;

    assert(false == ncbi_entrez.entrez_fetch(entrez_input, entrez_results));
    FS_dprint("\tUT_NCBIEntrez_00 - CASE_04 PASS");
    std::cout << "\tUT_NCBIEntrez_00 - CASE_04 PASS" << std::endl;

    FS_dprint("UT_NCBIEntrez_00 COMPLETE");
    std::cout << "UT_NCBIEntrez_00 COMPLETE" << std::endl;
}

void UnitTests::TestNCBIDatabase() {

    UTNCBIDatabase_00();

}

void UnitTests::UTNCBIDatabase_00() {
    NCBIDatabase ncbi_database = NCBIDatabase(mpFileSystem);
    NCBIDataResults_t ncbi_data_results;
    vect_str_t ncbi_data_accessions;

    FS_dprint("UT_NCBIDatabase_00 START");
    std::cout << "UT_NCBIDatabase_00 START" << std::endl;

    // Test 1, pulling gene id from protein database for a single NCBI ID
    std::string test_ncbi_id = "XP_014245616";
    std::string test_expected_geneid = "106664428";
    ncbi_data_accessions = {test_ncbi_id};
    ncbi_data_results = ncbi_database.get_ncbi_data(ncbi_data_accessions);

    assert(!ncbi_data_results.empty());
    assert(ncbi_data_results.find(test_ncbi_id) != ncbi_data_results.end());
    assert(ncbi_data_results.at(test_ncbi_id).geneid == test_expected_geneid);
    FS_dprint("\tUT_NCBIDatabase_00 - CASE_01 PASS");
    std::cout << "\tUT_NCBIDatabase_00 - CASE_01 PASS" << std::endl;


    // Test 2, pulling gene id from protein database for multiple NCBI ID
    //  including NCBI versions (XXXX.1)
    ncbi_data_accessions = {"XP_014245616", "XP_020482136.1", "XP_026476231.1"};
    ncbi_data_results = {};
    ncbi_data_results = ncbi_database.get_ncbi_data(ncbi_data_accessions);

    assert(!ncbi_data_results.empty());
    assert(ncbi_data_results.find("XP_014245616") != ncbi_data_results.end());
    assert(ncbi_data_results.at("XP_014245616").geneid == "106664428");
    assert(ncbi_data_results.find("XP_020482136.1") != ncbi_data_results.end());
    assert(ncbi_data_results.at("XP_020482136.1").geneid == "109976342");
    assert(ncbi_data_results.find("XP_026476231.1") != ncbi_data_results.end());
    assert(ncbi_data_results.at("XP_026476231.1").geneid == "113381705");
    FS_dprint("\tUT_NCBIDatabase_00 - CASE_02 PASS");
    std::cout << "\tUT_NCBIDatabase_00 - CASE_02 PASS" << std::endl;

    // Test 3, NCBIDatabase returns empty data if cannot find ID
    ncbi_data_accessions = {"XP_014245616safdsafs"};
    ncbi_data_results = {};

    assert(ncbi_database.get_ncbi_data(ncbi_data_accessions).empty());
    FS_dprint("\tUT_NCBIDatabase_00 - CASE_03 PASS");
    std::cout << "\tUT_NCBIDatabase_00 - CASE_03 PASS" << std::endl;

    // Test 4, NCBIDatabase returns empty data with no UID list
    ncbi_data_results = {};
    ncbi_data_accessions = {};

    assert(ncbi_database.get_ncbi_data(ncbi_data_accessions).empty());
    FS_dprint("\tUT_NCBIDatabase_00 - CASE_04 PASS");
    std::cout << "\tUT_NCBIDatabase_00 - CASE_04 PASS" << std::endl;

    FS_dprint("UT_NCBIDatabase_00 COMPLETE");
    std::cout << "UT_NCBIDatabase_00 COMPLETE" << std::endl;


}


#endif
