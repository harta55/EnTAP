/*
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

#include <csv.h>
#include "EntapDatabase.h"

/**
 * ======================================================================
 * Function EntapDatabase::EntapDatabase(FileSystem* filesystem)
 *
 * Description          - EnTAP Database constructor
 *                      - Controls all database related operations
 *                              - Downloading
 *                              - Indexing/Serializing
 *                              - Accession
 *                              - Creation
 *
 * Notes                - None
 *
 * @param filesystem    - Pointer to filesystem handler
 * @return              - None
 * ======================================================================
 */
EntapDatabase::EntapDatabase(FileSystem* filesystem) {
    // Initialize
    _pFilesystem     = filesystem;
    _temp_directory  = filesystem->get_temp_outdir();    // created previously
    _pSerializedDatabase = nullptr;
    _pDatabaseHelper     = nullptr;
    _use_serial          = true;                         // default
    _err_msg             = "";
    _err_code            = ERR_DATA_OK;
}


/**
 * ======================================================================
 * Function bool EntapDatabase::set_database(DATABASE_TYPE type)
 *
 * Description          - Initializes the chosen database (SQL/Serial)
 *
 * Notes                - None
 *
 * @param type          - Type of database to initialize
 * @param path          - Path to database
 * @return              - True/False is successful or not
 * ======================================================================
 */
bool EntapDatabase::set_database(DATABASE_TYPE type) {

    switch (type) {
        case ENTAP_SERIALIZED:
            // Filepath checked in routine
            _use_serial = true;
            return serialize_database_read(SERIALIZE_DEFAULT, ENTAP_DATABASE_BIN_PATH) == ERR_DATA_OK;
        case ENTAP_SQL:
            _use_serial = false;
            if (!_pFilesystem->file_exists(ENTAP_DATABASE_SQL_PATH)) {
                set_err_msg("Database not found at: "+ ENTAP_DATABASE_SQL_PATH, ERR_DATA_SET);
                return false;
            }
            if (_pDatabaseHelper != nullptr) return true;   // already generated
            _pDatabaseHelper = new SQLDatabaseHelper();
            return _pDatabaseHelper->open(ENTAP_DATABASE_SQL_PATH);
        default:
            return false;
    }
}


/**
 * ======================================================================
 * Function EntapDatabase::DATABASE_ERR EntapDatabase::download_database(
 *                                           EntapDatabase::DATABASE_TYPE type,
 *                                           std::string &path)
 *
 * Description          - Downloades the selected database from TreeGenes FTP
 *
 * Notes                - None
 *
 * @param type          - Type of database to download
 * @param path          - Path to output database
 * @return              - DATABASE_ERR type
 * ======================================================================
 */
EntapDatabase::DATABASE_ERR EntapDatabase::download_database(EntapDatabase::DATABASE_TYPE type, std::string &path) {
    DATABASE_ERR err;

    switch (type) {
        case ENTAP_SQL:
            err = download_entap_sql(path);
            break;
        case ENTAP_SERIALIZED:
            err = download_entap_serial(path);
            break;
        default:
            return ERR_DATA_OK;
    }
    _pFilesystem->delete_dir(_temp_directory);
    _pFilesystem->create_dir(_temp_directory);
    if (err != ERR_DATA_OK) _pFilesystem->delete_file(path);
    _err_code = err;
    return err;
}


/**
 * ======================================================================
 * Function EntapDatabase::DATABASE_ERR EntapDatabase::generate_database(
                            EntapDatabase::DATABASE_TYPE type,
                            std::string &path)
 *
 * Description          - Generates the EnTAP database (sql/serial) from online sources
 *
 * Notes                - None
 *
 * @param type          - Type of database to download
 * @param path          - Path to output database
 * @return              - DATABASE_ERR type
 * ======================================================================
 */
EntapDatabase::DATABASE_ERR EntapDatabase::generate_database(
                            EntapDatabase::DATABASE_TYPE type,
                            std::string &path) {
    DATABASE_ERR err;

    err = generate_entap_database(type, path);

    _pFilesystem->delete_dir(_temp_directory);
    _pFilesystem->create_dir(_temp_directory);
    if (err != (ERR_DATA_OK)) _pFilesystem->delete_file(path);
    _err_code = err;
    return err;
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_database(DATABASE_TYPE type, std::string &outpath) {
    DATABASE_ERR err_code;

    FS_dprint("Database type: " + ENTAP_DATABASE_TYPES_STR[type] + ", Outpath: " + outpath);
    switch (type) {

        case ENTAP_SQL:
            FS_dprint("Generating EnTAP SQL Database...");
            _use_serial = false;
            if (_pDatabaseHelper != nullptr) {
                // SQL should not already have been created
                return ERR_DATA_SQL_DUPLICATE;
            } else if (_pFilesystem->file_exists(outpath)) {
                // Out file should NOT exist yet
                set_err_msg("EnTAP SQL Database already exists at: " + outpath, ERR_DATA_FILE_EXISTS);
                return ERR_DATA_OK;
            }

            // Create sql database (will create if doesn't exist)
            _pDatabaseHelper = new SQLDatabaseHelper();
            if (!_pDatabaseHelper->create(outpath)) {
                // Return error if can't create
                set_err_msg("Unable to create the EnTAP SQL database", ERR_DATA_SQL_CREATE_DATABASE);
                return ERR_DATA_SQL_CREATE_DATABASE;
            }
            FS_dprint("Success!");
            break;

        case ENTAP_SERIALIZED:
            FS_dprint("Generating EnTAP Serialized Database...");
            _use_serial = true;
            if (_pSerializedDatabase != nullptr) {
                // Database should not already have been created
                set_err_msg("Serialized database already set!", ERR_DATA_SERIAL_DUPLICATE);
                return ERR_DATA_SERIAL_DUPLICATE;
            } else if (_pFilesystem->file_exists(outpath)) {
                // Out file should NOT exist yet
                set_err_msg("Serialized EnTAP database already exists at: " + outpath, ERR_DATA_FILE_EXISTS);
                return ERR_DATA_OK;
            }

            // Create serialized database
            _pSerializedDatabase = new EntapDatabaseStruct();
            FS_dprint("Success!");
            break;

        default:
            set_err_msg("ERROR: Unknown database type", ERR_DATA_UNHANDLED_TYPE);
            return ERR_DATA_UNHANDLED_TYPE;
    }

    // ---------------------- Add Database Entries ---------------------- //
    FS_dprint("Adding entries to database...");
    // Generate tax entries, don't need a path - using SQL member
    err_code = generate_entap_tax(type);
    if (err_code != ERR_DATA_OK) {
        return err_code;
    }

    // Generate go entries
    err_code = generate_entap_go(type);
    if (err_code != ERR_DATA_OK) {
        return err_code;
    }

    // Generate UniProt entries (this references GO database, must be done after)
    err_code = generate_entap_uniprot(type);
    if (err_code != ERR_DATA_OK) {
        return err_code;
    }
    // ------------------------------------------------------------------ //

    // Write database to file if necessary and set version number
    FS_dprint("All entries have been added, finalizing...");

    set_database_versions(type);
    switch (type) {
        case ENTAP_SQL:
            break;

        case ENTAP_SERIALIZED:
            FS_dprint("All entries added to database, serializing...");
            err_code = serialize_database_save(SERIALIZE_DEFAULT, outpath);
            if (err_code != ERR_DATA_OK) {
                FS_dprint("Unable to serialize database!");
                return err_code;
            }
            break;

        default:
            return ERR_DATA_MEM_ALLOC;
    }
    FS_dprint("Success! database generated to: " + outpath);
    return ERR_DATA_OK;
}

EntapDatabase::~EntapDatabase() {
    FS_dprint("Killing Object - EntapDatabase");
    if (_pDatabaseHelper != nullptr) {
        _pDatabaseHelper->close();
        delete(_pDatabaseHelper);
    }
    delete _pSerializedDatabase;
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_tax(EntapDatabase::DATABASE_TYPE type) {
    std::string temp_outpath;   // Path to unzipped files
    std::string sql_cmd;
    std::stringstream ss_temp;  // just for now
    std::string line;
    std::string ncbi_names_path;    // Path to uncompressed names file
    std::string ncbi_nodes_path;
    vect_str_t  split_line;         // Split line by tabs

    std::string tax_id;
    std::string tax_name;
    std::string lineage;

    // logging counts
    uint64 total_entries=0;
    uint64 current_entries=0;
    uint16 percent_complete;
    uint16 percent_prev=0;

    // Instead of creating tree, just using map to access each node
    // Keyed to NCBI ID's (first column of uncompressed files)
    std::unordered_map<std::string, TaxonomyNode> taxonomy_nodes;

    FS_dprint("Generating EnTAP Tax database entries...");

    if (type == ENTAP_SQL) {
        // If SQL, we want to create taxonomy table in database
        if (!create_sql_table(ENTAP_TAXONOMY)) {
            set_err_msg("Error generating SQL Taxonomy table", ERR_DATA_SQL_TAX_CREATE_TABLE);
            return ERR_DATA_SQL_TAX_CREATE_TABLE;
        }
    }

    temp_outpath = PATHS(_temp_directory, NCBI_TAX_DUMP_FILENAME);

    // download files from NCBI tax
    if (!_pFilesystem->download_ftp_file(FTP_NCBI_TAX_DUMP_TARGZ, temp_outpath)) {
        set_err_msg("Unable to download NCBI Taxonomy FTP files" + _pFilesystem->get_error(), ERR_DATA_TAX_DOWNLOAD);
        return ERR_DATA_TAX_DOWNLOAD;
    }
    // decompress TAR.GZ file
    if (!_pFilesystem->decompress_file(temp_outpath, _temp_directory, FileSystem::FILE_TAR_GZ)) {
        set_err_msg("Unable to decompress NCBI Taxonomy data" + _pFilesystem->get_error(), ERR_DATA_FILE_DECOMPRESS);
        return ERR_DATA_FILE_DECOMPRESS;
    }
    // remove tar file
    _pFilesystem->delete_file(temp_outpath);
    // set paths to uncompressed files
    ncbi_names_path = PATHS(_temp_directory, NCBI_TAX_DUMP_FTP_NAMES);
    ncbi_nodes_path = PATHS(_temp_directory, NCBI_TAX_DUMP_FTP_NODES);

    FS_dprint("Files downloaded and compressed, parsing...");


    FS_dprint("Parsing NCBI Names file at: " + ncbi_names_path);
    try {
        // parse through names of taxonomy ID's and add to map
        std::ifstream infile(ncbi_names_path);
        while (std::getline(infile, line)) {
            if (line == "") continue;
            total_entries++;
            // Split line by tabs, lazy..fix :(
            split_line = split_string(line, NCBI_TAX_DUMP_DELIM);

            tax_id = split_line[NCBI_TAX_DUMP_COL_ID];
            tax_name   = split_line[NCBI_TAX_DUMP_COL_NAME];

            // Check if map already has this entry, if not - generate
            if (taxonomy_nodes.find(tax_id) == taxonomy_nodes.end()) {
                taxonomy_nodes.emplace(tax_id, TaxonomyNode(tax_id));
            }

            std::unordered_map<std::string, TaxonomyNode>::iterator it = taxonomy_nodes.find(tax_id);

            // We'll want to use scientific names when displaying lineage
            if (split_line[NCBI_TAX_DUMP_COL_NAME_CLASS].compare(NCBI_TAX_DUMP_SCIENTIFIC) ==0) {
                it->second.sci_name = tax_name;
            }
            it->second.names.push_back(tax_name);
        }
        infile.close();
        FS_dprint("Success! Parsing nodes file at: " + ncbi_nodes_path);

        // parse through nodes file
        std::ifstream infile_node(ncbi_nodes_path);
        while (std::getline(infile_node, line)) {
            if (line == "") continue;

            // split line by tabs
            split_line = split_string(line, NCBI_TAX_DUMP_DELIM);

            tax_id = split_line[NCBI_TAX_DUMP_COL_ID];

            // Ensure node has name associated with it
            std::unordered_map<std::string, TaxonomyNode>::iterator it = taxonomy_nodes.find(tax_id);

            // Node with no names, skip we don't want this
            if (it == taxonomy_nodes.end()) continue;

            // Set parent node NCBI ID
            it->second.parent_id = split_line[NCBI_TAX_DUMP_COL_PARENT];
        }
        infile_node.close();
        FS_dprint("Success! Compiling final NCBI results...");

        // parse through entire map and generate NCBI taxonomy entries
        TaxEntry taxEntry;
        for (auto &pair : taxonomy_nodes) {
            // want a separate entry for each name (doing this for now, may change)
            for (std::string name : pair.second.names) {
                LOWERCASE(name);
                current_entries++;
                taxEntry = {};
                taxEntry.lineage = entap_tax_get_lineage(pair.second, taxonomy_nodes);
                LOWERCASE(taxEntry.lineage);
                taxEntry.tax_id  = pair.second.ncbi_id;
                taxEntry.tax_name= name;

                // Add to SQL database or other...
                if (type == ENTAP_SQL) {
                    if (!sql_add_tax_entry(taxEntry)) {
                        // unable to add entry
                        set_err_msg("Unable to add Taxonomy entry " + name, ERR_DATA_SQL_CREATE_ENTRY);
                        return ERR_DATA_SQL_CREATE_ENTRY;
                    }
                } else {
                    // Add to database that will be serialized
                    _pSerializedDatabase->taxonomic_data[name] = taxEntry;
                }
            }
            // ********************* logging **************************** //
            percent_complete = (uint16) round((fp32)current_entries / total_entries * 100);
            if (percent_complete % STATUS_UPDATES == 0 && percent_complete != percent_prev) {
                FS_dprint("Percent complete: "+ std::to_string(percent_complete) + "%");
                percent_prev = percent_complete;
            }
            // ********************************************************** //
        }
    } catch (const std::exception &e) {
        set_err_msg("Unable to parse taxonomy data: " + std::string(e.what()), ERR_DATA_TAXONOMY_PARSE);
        return ERR_DATA_TAXONOMY_PARSE;
    }

    FS_dprint("Success! NCBI data complete");
    return ERR_DATA_OK;
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_go(EntapDatabase::DATABASE_TYPE type) {
    FS_dprint("Generating EnTAP Gene Ontology entries...");

    std::string go_db_path;
    std::string go_term_path;
    std::string go_graph_path;
    std::string go_database_targz;  // Outpath to downloaded tar.gz file
    std::string go_database_dir;    // Directory that will contain go files

    go_database_targz = PATHS(_temp_directory, GO_TERMDB_FILE);


    // download Gene Ontology database file
    if (!_pFilesystem->download_ftp_file(FTP_GO_DATABASE, go_database_targz)) {
        // failed to download from ftp
        set_err_msg("Unable to download GO data from FTP address " + FTP_GO_DATABASE, ERR_DATA_GO_DOWNLOAD);
        return ERR_DATA_GO_DOWNLOAD;
    }

    // decompress database file
    if (!_pFilesystem->decompress_file(go_database_targz, _temp_directory, FileSystem::FILE_TAR_GZ)) {
        // failed to decompress
        set_err_msg("Unable to decompress GO database file " + _pFilesystem->get_error(), ERR_DATA_GO_DECOMPRESS);
        return ERR_DATA_GO_DECOMPRESS;
    }
    _pFilesystem->delete_file(go_database_targz);

    // Files are packaged within a directory. Set paths
    go_database_dir = PATHS(_temp_directory, GO_TERMDB_DIR);
    go_term_path    = PATHS(go_database_dir, GO_TERM_FILE);
    go_graph_path   = PATHS(go_database_dir, GO_GRAPH_FILE);

    if (!_pFilesystem->file_exists(go_term_path) || !_pFilesystem->file_exists(go_graph_path)) {
        set_err_msg("Necessary Gene Ontology files do not exist at:\n" + go_term_path +
                   "\n" + go_graph_path, ERR_DATA_GO_DOWNLOAD);
        return ERR_DATA_GO_DOWNLOAD;
    }

    // If we are creating SQL database, add GO table
    if (type == ENTAP_SQL) {
        if (!create_sql_table(ENTAP_GENE_ONTOLOGY)) {
            // error creating table
            return ERR_DATA_SQL_GO_CREATE_TABLE;
        }
    }

    try {
        // Parse through graph file
        io::CSVReader<6, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(go_graph_path);
        std::string index,root,branch, temp, distance, temp2;
        std::map<std::string,std::string> distance_map;
        while (in.read_row(index,root,branch, temp, distance, temp2)) {
            if (root.compare(GO_BIOLOGICAL_LVL) == 0     ||
                root.compare(GO_MOLECULAR_LVL) == 0  ||
                root.compare(GO_CELLULAR_LVL) ==0) {
                if (distance_map.find(branch) == distance_map.end()) {
                    distance_map.emplace(branch,distance);
                } else {
                    if (distance.empty()) continue;
                    fp32 cur   = std::stoi(distance_map[branch]);
                    fp32 query = std::stoi(distance);
                    if (query > cur) distance_map[branch] = distance;
                }
            }
        }
        GoEntry goEntry;
        std::string num,term,cat,go,ex,ex1,ex2;
        io::CSVReader<7, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in2(go_term_path);
        while (in2.read_row(num,term,cat,go,ex,ex1,ex2)) {
            goEntry = {};
            goEntry.category = cat;
            goEntry.level = distance_map[num];
            goEntry.term = term;
            goEntry.go_id = go;

            // Add to SQL database OR to overall map
            if (type == ENTAP_SQL) {
                if (!sql_add_go_entry(goEntry)) {
                    set_err_msg("Unable to add GO entry: " + goEntry.go_id, ERR_DATA_GO_ENTRY);
                    return ERR_DATA_GO_ENTRY;
                }
            } else {
                _pSerializedDatabase->gene_ontology_data[go] = goEntry;
            }
        }
    } catch (const std::exception &e) {
        set_err_msg("Unable to parse Gene Ontology data: " + std::string(e.what()), ERR_DATA_GO_PARSE);
        return ERR_DATA_GO_PARSE;
    }

    FS_dprint("Success! Gene Ontology data complete");
    _pFilesystem->delete_file(go_graph_path);
    _pFilesystem->delete_file(go_term_path);
    return ERR_DATA_OK;
}


EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_uniprot(EntapDatabase::DATABASE_TYPE type) {

    std::string uniprot_flat_gz;
    std::string uniprot_flat;       // Decompressed path (this is parsed)
    std::string line;
    std::string dat_tag;
    std::string database;
    std::string go_list;            // Will be turned into go_format when indexed (comma delim)
    std::string kegg_list;
    uint16      file_status;
    uint16      index_go;
    uint16      index_term;
    std::string data;
    bool        same_entry = false;

    // Set output path for FTP file
    uniprot_flat_gz = PATHS(_temp_directory, UNIPROT_DAT_FILE_GZ);
    uniprot_flat    = PATHS(_temp_directory, UNIPROT_DAT_FILE);


    // download UniProt flat file
    if (!_pFilesystem->download_ftp_file(FTP_UNIPROT_FLAT_FILE, uniprot_flat_gz)) {
        // failed to download from ftp
        set_err_msg("Unable to download UniProt data from " + FTP_UNIPROT_FLAT_FILE + _pFilesystem->get_error(), ERR_DATA_UNIPROT_DOWNLOAD);
        return ERR_DATA_UNIPROT_DOWNLOAD;
    }

    // decompress database file
    if (!_pFilesystem->decompress_file(uniprot_flat_gz, uniprot_flat, FileSystem::FILE_GZ)) {
        // failed to decompress
        set_err_msg("Unable to decompress UniProt data" + _pFilesystem->get_error(), ERR_DATA_UNIPROT_DECOMPRESS);
        return ERR_DATA_UNIPROT_DECOMPRESS;
    }
    _pFilesystem->delete_file(uniprot_flat_gz); // Ensure gz file is deleted

    FS_dprint("UniProt file downloaded/decompressed, verifying...");
    // Ensure file is valid (should be)
    file_status = _pFilesystem->get_file_status(uniprot_flat);
    if (file_status != 0) {
        // File invalid
        set_err_msg(_pFilesystem->print_file_status(file_status, uniprot_flat), ERR_DATA_UNIPROT_FILE);
        return ERR_DATA_UNIPROT_FILE;
    }

    // If we are creating SQL database, add UniProt table
    if (type == ENTAP_SQL) {
        if (!create_sql_table(ENTAP_UNIPROT)) {
            // error creating table
            set_err_msg("Unable to create UniProt SQL Table", ERR_DATA_SQL_UNIPROT_CREATE_TABLE);
            return ERR_DATA_SQL_UNIPROT_CREATE_TABLE;
        }
    }
    FS_dprint("UniProt file successfully downloaded/decompressed. Parsing...");
    // File valid, continue to parse
    try {
        UniprotEntry uniprotEntry;
        std::ifstream infile(uniprot_flat);
        while (std::getline(infile, line)) {
            if (line.empty() || line.length() < UNIPROT_DAT_TAG_LEN) continue;
            STR_ERASE(line, '\n');
            dat_tag = line.substr(0, UNIPROT_DAT_TAG_LEN);
            if (line.length() > UNIPROT_DAT_TAG_DATA_POS) {
                data = line.substr(UNIPROT_DAT_TAG_DATA_POS); // All data from the line
            }

            // Check if this is sequence ID
            if (dat_tag == UNIPROT_DAT_TAG_ID) {
                if (same_entry) {
                    FS_dprint("ERROR: Same entry is true");
                    return ERR_DATA_UNIPROT_PARSE;
                }
                same_entry = true;
                uniprotEntry.uniprot_id = line.substr(UNIPROT_DAT_TAG_DATA_POS,
                                                      line.find_first_of(' ',UNIPROT_DAT_TAG_DATA_POS)-UNIPROT_DAT_TAG_DATA_POS);
            } else if (dat_tag == UNIPROT_DAT_TAG_DATABASE_X_REF) {
                // ID matches the database cross reference ID

                // Check which database we have
                database = data.substr(0, data.find(UNIPROT_DAT_TAG_DATABASE_DELIM));
                if (database == UNIPROT_DAT_TAG_DATABASE_GO) {
                    index_go = (uint16) data.find("GO:");
                    index_term = (uint16) data.find(';', index_go);
                    go_list += data.substr(index_go, index_term - index_go) + ',';

                } else if (database == UNIPROT_DAT_TAG_DATABASE_KEGG) {
                    index_go = (uint16) data.find("KEGG:");
                    index_term = (uint16) data.find(';', index_go);
                    kegg_list += data.substr(index_go, index_term - index_go) + ',';
                } else {
                    // Neither GO nor KEGG, add to x refs
                    uniprotEntry.database_x_refs += "|" + data;
                }

            } else if (dat_tag == UNIPROT_DAT_TAG_COMMENT) {
                // ID matches the comments (sometimes useful info? maybe)
                uniprotEntry.comments += "|" + data;
            } else if (dat_tag == UNIPROT_DAT_TAG_NEXT_ENTRY){
                // We've hit the next entry, add previous to the database
                if (go_list.length() > 0) go_list.pop_back();   // remove trailing ','
                uniprotEntry.go_terms = format_go_delim(go_list, ',');
                uniprotEntry.kegg_terms = kegg_list;

                if (!add_uniprot_entry(type, uniprotEntry)) {
                    // Unable to add entry
                    set_err_msg("ERROR: Unable to add entry:\n" + uniprotEntry.print(), ERR_DATA_UNIPROT_ENTRY);
                    return ERR_DATA_UNIPROT_ENTRY;
                }
                uniprotEntry = {};
                same_entry = false;
                go_list = "";
                kegg_list = "";

            } else {
                // Unhandled information from UniProt mapping, discard
            }
        }
    } catch (const std::exception &e) {
        set_err_msg("Unable to parse UniProt data: " + std::string(e.what()), ERR_DATA_UNIPROT_PARSE);
        return ERR_DATA_UNIPROT_PARSE;
    }

    FS_dprint("Success! UniProt entries added");
    _pFilesystem->delete_file(uniprot_flat);
    return ERR_DATA_OK;
}


// WARNING: recursive
std::string EntapDatabase::entap_tax_get_lineage(EntapDatabase::TaxonomyNode &node,
                                                 std::unordered_map<std::string,TaxonomyNode>& map) {
    if (node.ncbi_id == "1" || node.ncbi_id == "") {
        return "root";
    } else {
        return node.sci_name + ";" + entap_tax_get_lineage(map.at(node.parent_id), map);
    }
}

bool EntapDatabase::sql_add_tax_entry(TaxEntry &taxEntry) {
    char *sql_cmd;

    if (_pDatabaseHelper == nullptr) return false;

    sql_cmd = sqlite3_mprintf(
            "INSERT INTO %Q (%Q,%Q,%Q) "\
            "VALUES (%Q, %Q, %Q);",

            SQL_TABLE_NCBI_TAX_TITLE.c_str(),
            SQL_COL_NCBI_TAX_TAXID.c_str(),
            SQL_COL_NCBI_TAX_LINEAGE.c_str(),
            SQL_COL_NCBI_TAX_NAME.c_str(),
            taxEntry.tax_id.c_str(),
            taxEntry.lineage.c_str(),
            taxEntry.tax_name.c_str()
    );
    return _pDatabaseHelper->execute_cmd(sql_cmd);
}

bool EntapDatabase::sql_add_go_entry(GoEntry &goEntry) {
    char *sql_cmd;

    if (_pDatabaseHelper == nullptr) return false;

    sql_cmd = sqlite3_mprintf(
            "INSERT INTO %Q (%Q,%Q,%Q,%Q) "\
            "VALUES (%Q, %Q, %Q, %Q);",

            SQL_TABLE_GO_TITLE.c_str(),
            SQL_TABLE_GO_COL_ID.c_str(),
            SQL_TABLE_GO_COL_DESC.c_str(),
            SQL_TABLE_GO_COL_CATEGORY.c_str(),
            SQL_TABLE_GO_COL_LEVEL.c_str(),
            goEntry.go_id.c_str(),
            goEntry.term.c_str(),
            goEntry.category.c_str(),
            goEntry.level.c_str()
    );
    return _pDatabaseHelper->execute_cmd(sql_cmd);
}

bool EntapDatabase::add_uniprot_entry(EntapDatabase::DATABASE_TYPE type, UniprotEntry &entry) {
    bool ret = true;
    char *sql_cmd;

    switch (type) {
        case ENTAP_SERIALIZED:
            if (_pSerializedDatabase == nullptr) {
                FS_dprint("ERROR: Serialized database NULL");
                ret = false;
                break;
            }
            _pSerializedDatabase->uniprot_data[entry.uniprot_id] = entry;
            break;

        case ENTAP_SQL:
            if (_pDatabaseHelper == nullptr) {
                FS_dprint("ERROR: SQL database NULL");
                ret = false;
                break;
            }

            sql_cmd = sqlite3_mprintf(
                    "INSERT INTO %Q (%Q,%Q,%Q) " \
                    "VALUES (%Q,%Q,%Q);",

                    SQL_TABLE_UNIPROT_TITLE.c_str(),
                    SQL_TABLE_UNIPROT_COL_ID.c_str(),
                    SQL_TABLE_UNIPROT_COL_XREF.c_str(),
                    SQL_TABLE_UNIPROT_COL_COMM.c_str(),
                    entry.uniprot_id.c_str(),
                    entry.database_x_refs.c_str(),
                    entry.comments.c_str()
            );
            ret = _pDatabaseHelper->execute_cmd(sql_cmd);
            break;

        default:
            break;
    }
    return ret;
}

bool EntapDatabase::create_sql_table(DATABASE_TYPE type) {
    char *sql_cmd;
    bool success;

    if (_pDatabaseHelper == nullptr) return false;

    switch (type) {
        case ENTAP_TAXONOMY:
            FS_dprint("Creating SQL Taxonomy table...");
            sql_cmd = sqlite3_mprintf(
                    "CREATE TABLE %Q ("                 \
                    "ID      INTEGER PRIMARY KEY     NOT NULL," \
                    "%Q      TEXT                NOT NULL," \
                    "%Q      TEXT                NOT NULL," \
                    "%Q      TEXT                NOT NULL);",

                    SQL_TABLE_NCBI_TAX_TITLE.c_str(),
                    SQL_COL_NCBI_TAX_TAXID.c_str(),
                    SQL_COL_NCBI_TAX_LINEAGE.c_str(),
                    SQL_COL_NCBI_TAX_NAME.c_str()
            );
            break;

        case ENTAP_GENE_ONTOLOGY:
            FS_dprint("Creating SQL Gene Ontology table...");
            sql_cmd = sqlite3_mprintf(
                    "CREATE TABLE %Q ("                     \
                    "ID        INTEGER PRIMARY KEY       NOT NULL," \
                    "%Q        TEXT                      NOT NULL," \
                    "%Q        TEXT                      NOT NULL," \
                    "%Q        TEXT                      NOT NULL," \
                    "%Q        TEXT                      NOT NULL);",

                    SQL_TABLE_GO_TITLE.c_str(),
                    SQL_TABLE_GO_COL_ID.c_str(),
                    SQL_TABLE_GO_COL_DESC.c_str(),
                    SQL_TABLE_GO_COL_CATEGORY.c_str(),
                    SQL_TABLE_GO_COL_LEVEL.c_str()
            );
            break;

        case ENTAP_UNIPROT:
            FS_dprint("Creating SQL UniProt table...");
            sql_cmd = sqlite3_mprintf(
                    "CREATE TABLE %Q ("                 \
                    "ID    INTEGER PRIMARY KEY     NOT NULL," \
                    "%Q    TEXT                    NOT NULL," \
                    "%Q    TEXT                    NOT NULL," \
                    "%Q    TEXT                    NOT NULL);",

                    SQL_TABLE_UNIPROT_TITLE.c_str(),
                    SQL_TABLE_UNIPROT_COL_ID.c_str(),
                    SQL_TABLE_UNIPROT_COL_XREF.c_str(),
                    SQL_TABLE_UNIPROT_COL_COMM.c_str()
            );
            break;

        case ENTAP_VERSION:
            FS_dprint("Creating version table...");
            sql_cmd = sqlite3_mprintf(
                    "CREATE TABLE %Q ("                     \
                    "ID   INTEGER PRIMARY KEY    NOT NULL," \
                    "%Q   TEXT                   NOT NULL);",
                    SQL_TABLE_VERSION_TITLE.c_str(),
                    SQL_TABLE_VERSION_COL_VER.c_str()
                    );
            break;

        default:
            FS_dprint("ERROR: Unhandled SQL table creation");
            return false;
    }

    success = _pDatabaseHelper->execute_cmd(sql_cmd);
    if (success) {
        FS_dprint("Success!");
    } else {
        std::string temp = sql_cmd;
        FS_dprint("ERROR: Unable to create table with command: \n" + temp);
    }
    return success;
}

EntapDatabase::DATABASE_ERR EntapDatabase::download_entap_serial(std::string &out_path) {
    std::string temp_gz_path;

    FS_dprint("Downloading EnTAP serialized database...");

    // set temp path (will be downloaded to this then decompressed)
    temp_gz_path = PATHS(_temp_directory, Defaults::ENTAP_DATABASE_SERIAL_GZ);

    // download file (will be compressed as gz)
    if (!_pFilesystem->download_ftp_file(FTP_ENTAP_DATABASE_SERIAL, temp_gz_path)) {
        // File download failed!
        set_err_msg("Unable to download EnTAP Serial Database from " + FTP_ENTAP_DATABASE_SERIAL +
            _pFilesystem->get_error(), ERR_DATA_SERIAL_FTP);
        return ERR_DATA_SERIAL_FTP;
    }

    // decompress file to outpath
    if (!_pFilesystem->decompress_file(temp_gz_path, out_path, FileSystem::FILE_GZ)) {
        // Decompression failed!
        set_err_msg("Unable to decompress EnTAP Serial Database at " + temp_gz_path + _pFilesystem->get_error(),
            ERR_DATA_SERIAL_DECOMPRESS);
        return ERR_DATA_SERIAL_DECOMPRESS;
    }

    // remove compressed file (already
    _pFilesystem->delete_file(temp_gz_path);

    return ERR_DATA_OK;
}

EntapDatabase::DATABASE_ERR EntapDatabase::download_entap_sql(std::string &path) {
    std::string temp_gz_path;

    FS_dprint("Downloading EnTAP sql database...");

    // set temp path (will be downloaded to this then decompressed)
    temp_gz_path = PATHS(_temp_directory, Defaults::ENTAP_DATABASE_SQL_GZ);

    // download file (will be compressed as gz)
    if (!_pFilesystem->download_ftp_file(FTP_ENTAP_DATABASE_SQL, temp_gz_path)) {
        // File download failed!
        set_err_msg("Unable to download EnTAP Serial Database from " + FTP_ENTAP_DATABASE_SQL +
                    _pFilesystem->get_error(), ERR_DATA_SQL_FTP);
        return ERR_DATA_SQL_FTP;
    }

    // decompress file to outpath
    if (!_pFilesystem->decompress_file(temp_gz_path, path, FileSystem::FILE_GZ)) {
        // Decompression failed!
        set_err_msg("Unable to decompress EnTAP Serial Database at " + temp_gz_path + _pFilesystem->get_error(),
                    ERR_DATA_SQL_DECOMPRESS);
        return ERR_DATA_SQL_DECOMPRESS;
    }

    // remove compressed file (already done)
    _pFilesystem->delete_file(temp_gz_path);

    return ERR_DATA_OK;
}

GoEntry EntapDatabase::get_go_entry(std::string &go_id) {
    GoEntry goEntry;

    if (go_id.empty()) return GoEntry();

    if (_use_serial) {
        // Using serialized database
        go_serial_map_t::iterator it = _pSerializedDatabase->gene_ontology_data.find(go_id);
        if (it == _pSerializedDatabase->gene_ontology_data.end()) {
//            FS_dprint("Unable to find GO ID: " + go_id);
            return GoEntry();
        } else return it->second;

    } else {
        // Using SQL database
        std::vector<std::vector<std::string>> results;
        // Check temp if previously found (increase speeds)
        go_serial_map_t::iterator it = _sql_go_helper.find(go_id);
        if (it != _sql_go_helper.end()) return it->second;
        // Generate SQL query
        char *query = sqlite3_mprintf(
                "SELECT %q, %q, %q, %q FROM %q WHERE %q=%Q",
                SQL_TABLE_GO_COL_ID.c_str(),
                SQL_TABLE_GO_COL_DESC.c_str(),
                SQL_TABLE_GO_COL_CATEGORY.c_str(),
                SQL_TABLE_GO_COL_LEVEL.c_str(),
                SQL_TABLE_GO_TITLE.c_str(),
                SQL_TABLE_GO_COL_ID.c_str(),
                go_id.c_str()
        );
        try {
            if (results.empty()) return GoEntry();
            results = _pDatabaseHelper->query(query);
            goEntry.go_id    = results[0][0];
            goEntry.term     = results[0][1];
            goEntry.category = results[0][2];
            goEntry.level    = results[0][3];
            _sql_go_helper[go_id] = goEntry;
            return goEntry;
        } catch (std::exception &e) {
            // Do not fatal error
            FS_dprint(e.what());
            return GoEntry();
        }
    }
}

TaxEntry EntapDatabase::get_tax_entry(std::string &species) {
    TaxEntry taxEntry;
    std::string temp_species;
    uint64 index;

    if (species.empty()) return TaxEntry();

    LOWERCASE(species); // ensure lowercase (database is based on this for direct matching)

    if (_use_serial) {
        // Using serialized database
        tax_serial_map_t::iterator it = _pSerializedDatabase->taxonomic_data.find(species);
        if (it == _pSerializedDatabase->taxonomic_data.end()) {
            // If we can't find species, keep trying by making it more broad
            temp_species = species;
            while (true) {
                index = temp_species.find_last_of(" ");
                if (index == std::string::npos) break;
                temp_species = temp_species.substr(0, index);
                it = _pSerializedDatabase->taxonomic_data.find(temp_species);
                if (it != _pSerializedDatabase->taxonomic_data.end()) {
                    return it->second;
                }
            }
            return TaxEntry();
        } else return it->second;

    } else {
        // Using SQL database
        std::vector<std::vector<std::string>> results;
        temp_species = species;
        try {
            // If we can't find species, keep trying by making it more broad
            while (true) {
                // Generate SQL query
                char *query = sqlite3_mprintf(
                        "SELECT %q, %q FROM %q WHERE %q=%Q",
                        SQL_COL_NCBI_TAX_TAXID.c_str(),
                        SQL_COL_NCBI_TAX_LINEAGE.c_str(),
                        SQL_TABLE_NCBI_TAX_TITLE.c_str(),
                        SQL_COL_NCBI_TAX_NAME.c_str(),
                        temp_species.c_str()
                );
                results = _pDatabaseHelper->query(query);
                if (results.empty()) {
                    index = temp_species.find_last_of(" ");
                    if (index == std::string::npos) return TaxEntry(); // couldn't find
                    temp_species = temp_species.substr(0, index);
                } else break; // Found species
            }

            taxEntry.tax_id  = results[0][0];
            taxEntry.lineage = results[0][1];
            taxEntry.tax_name= temp_species;
            return taxEntry;

        } catch (std::exception &e) {
            // Do not fatal error
            FS_dprint(e.what());
            return TaxEntry();
        }
    }
}

UniprotEntry EntapDatabase::get_uniprot_entry(std::string& accession) {
    UniprotEntry uniprotEntry;

    if (accession.empty()) return UniprotEntry();

    try {
        if (_use_serial) {
            // Using serialized database
            uniprot_serial_map_t::iterator it =
                    _pSerializedDatabase->uniprot_data.find(accession);
            if (it == _pSerializedDatabase->uniprot_data.end()) {
                FS_dprint("Unable to find Uniprot Entry: " + accession);
                return UniprotEntry();
            } else return it->second;
        } else {
            // Using SQL database
            std::vector<std::vector<std::string>> results;
            char *query = sqlite3_mprintf(
                    "SELECT %q, %q, %q FROM %q WHERE %q=%Q",
                    SQL_TABLE_UNIPROT_COL_ID.c_str(),
                    SQL_TABLE_UNIPROT_COL_XREF.c_str(),
                    SQL_TABLE_UNIPROT_COL_COMM.c_str(),
                    SQL_TABLE_UNIPROT_TITLE.c_str(),
                    SQL_TABLE_UNIPROT_COL_ID.c_str(),
                    accession.c_str()
            );
            if (results.empty()) return UniprotEntry();
            results = _pDatabaseHelper->query(query);
            uniprotEntry.uniprot_id      = results[0][0];
            uniprotEntry.database_x_refs = results[0][1];
            uniprotEntry.comments        = results[0][2];
            return uniprotEntry;
        }
    } catch (const std::exception &e) {
        FS_dprint("ERROR: Unhandled finding Uniprot Entry: "+ std::string(e.what()));
        return UniprotEntry();
    }
}

EntapDatabase::DATABASE_ERR EntapDatabase::serialize_database_save(SERIALIZATION_TYPE type, std::string &out_path) {
    FS_dprint("Serializing EnTAP database to:" + out_path);

    if (_pSerializedDatabase == nullptr) {
        // Error in allocating memory
        FS_dprint("Error allocating memory to EnTAP Database");
        return ERR_DATA_SERIALIZE_SAVE;
    }

    try {
        FS_dprint("Serializing type: " + std::to_string(type));
        std::ofstream file(out_path);
        switch (type) {

#ifdef USE_BOOST
            case BOOST_TEXT_ARCHIVE: {
                boostAR::text_oarchive oa(file);
                oa << *_pSerializedDatabase;
                break;
            }

            case BOOST_BIN_ARCHIVE: {
                boostAR::binary_oarchive oa_bin(file);
                oa_bin << *_pSerializedDatabase;
                break;
            }
#else // Use CEREAL for serilization

            case CEREAL_BIN_ARCHIVE: {
                std::stringstream ss;
                cereal::BinaryOutputArchive oarchive(ss); // Create an output archive
                oarchive(*_pSerializedDatabase); // Write the data to the archive

                // Write string stream to file
                file << ss.str();
                break;
            }   // archive out of scope, flush
#endif

            default:
                return ERR_DATA_SERIALIZE_SAVE;
        }
        file.close();
    } catch (std::exception &e) {
        set_err_msg("Unable to serialize EnTAP database\n" + std::string(e.what()),
                    ERR_DATA_SERIALIZE_SAVE);
        return ERR_DATA_SERIALIZE_SAVE;
    }
    return ERR_DATA_OK;
}

EntapDatabase::DATABASE_ERR EntapDatabase::serialize_database_read(SERIALIZATION_TYPE type, std::string &in_path) {
    FS_dprint("Reading serialized database from: " + in_path + "\n Of type: " + std::to_string(type));

    if (!_pFilesystem->file_exists(in_path)) {
        set_err_msg("Serialized EnTAP database does not exist at: " + in_path, ERR_DATA_SERIALIZE_READ);
        return ERR_DATA_SERIALIZE_READ;
    }

    // Already generated? If no, generate
    if (_pSerializedDatabase != nullptr) return ERR_DATA_OK;
    _pSerializedDatabase = new EntapDatabaseStruct();

    try {
        switch (type) {
#ifdef USE_BOOST
            case BOOST_TEXT_ARCHIVE:
            {
                std::ifstream ifs(in_path);
                boost::archive::text_iarchive ia(ifs);
                ia >> *_pSerializedDatabase;
                ifs.close();
                break;
            }

            case BOOST_BIN_ARCHIVE:
            {
                std::ifstream ifs(in_path);
                boost::archive::binary_iarchive ia(ifs);
                ia >> *_pSerializedDatabase;
                ifs.close();
                break;
            }
#else
            case CEREAL_BIN_ARCHIVE: {
                std::stringstream ss;
                std::ifstream in_file(in_path);
                ss << in_file.rdbuf();
                in_file.close();
                cereal::BinaryInputArchive iarchive(ss); // Create an input archive
                iarchive(*_pSerializedDatabase); // Read the data from the archive
                break;
            }
#endif

            default:
                return ERR_DATA_SERIALIZE_READ;
        }

    } catch (const std::exception &e) {
        set_err_msg("Unable to read Serialized EnTAP database: " + std::string(e.what()), ERR_DATA_SERIALIZE_READ);
        return ERR_DATA_SERIALIZE_READ;
    }
    if (!is_valid_version()) {
        set_err_msg("EnTAP database version is not compatible with this version of EnTAP.\n" \
                    "Current Version: " + get_current_version_str() + "\nRequired version: " +
                            get_required_version_str(),
                    ERR_DATA_INCOMPATIBLE_VER);
        FS_dprint("WARNING: invalid database version!!!");
        return ERR_DATA_INCOMPATIBLE_VER;
    }
    return ERR_DATA_OK;
}

std::string EntapDatabase::print_error_log() {
    return "\nEnTAP Database Error: " + _err_msg;
}

void EntapDatabase::set_err_msg(std::string msg, DATABASE_ERR code) {
    FS_dprint(msg);
    _err_msg = msg;
    _err_code = code;
}


EntapDatabase::TaxonomyNode::TaxonomyNode(std::string id) {
    parent_id = "";
    sci_name = "";
    ncbi_id = id;
}

// terms = "GO:4321431,GO:807890"
go_format_t EntapDatabase::format_go_delim(std::string terms, char delim) {
    go_format_t output;
    std::string temp;
    std::vector<std::vector<std::string>>results;

    if (terms.empty()) return output;
    std::istringstream ss(terms);
    while (std::getline(ss,temp,delim)) {
        GoEntry term_info = get_go_entry(temp);
        if (!term_info.is_empty()) {
            output[term_info.category].push_back(temp + "-" + term_info.term +
                                                 "(L=" + term_info.level + ")");
        }
    }
    return output;
}

bool EntapDatabase::is_valid_version() {
    return (get_current_version_str() == get_required_version_str());
}

std::string EntapDatabase::get_current_version_str() {
    std::string version_str;

    if (_use_serial) {
        // Using serialized database
        if (_pSerializedDatabase != nullptr) {
            version_str = std::to_string(_pSerializedDatabase->MAJOR_VERSION) + "." +
                          std::to_string(_pSerializedDatabase->MINOR_VERSION);
        } else {
            version_str = "";
        }

    } else {
        // Using SQL database
        char* query;

        if (_pDatabaseHelper != nullptr) {
            query = sqlite3_mprintf(
                    "SELECT %Q FROM %Q",
                    SQL_TABLE_VERSION_TITLE.c_str(),
                    SQL_TABLE_VERSION_COL_VER.c_str()
            );

            try {
                FS_dprint("Executing sql cmd: " + std::string(query));
                version_str = _pDatabaseHelper->query(query)[0][0];
                FS_dprint("Success! Returned version: " + version_str);

            } catch (...) {
                FS_dprint("ERROR: couldn't get SQL version");
                version_str = "";
            }

        } else {
            version_str = "";
        }
    }

    return version_str;
}

std::string EntapDatabase::get_required_version_str() {
    if (_use_serial) {
        return std::to_string(SERIALIZE_MAJOR) + "." + std::to_string(SERIALIZE_MINOR);
    } else {
        return std::to_string(SQL_MAJOR) + "." + std::to_string(SQL_MINOR);
    }
}

bool EntapDatabase::set_database_versions(EntapDatabase::DATABASE_TYPE type) {
    bool ret;

    switch (type) {

        case ENTAP_SERIALIZED:
            if (_pSerializedDatabase != nullptr) {
                _pSerializedDatabase->MAJOR_VERSION = SERIALIZE_MAJOR;
                _pSerializedDatabase->MINOR_VERSION = SERIALIZE_MINOR;
                ret = true;
            } else {
                ret = false;
            }
            break;

        case ENTAP_SQL:
            if (_pDatabaseHelper != nullptr) {
                char *query;
                std::string version_str;

                // Can we create the table to store the versioning data
                if (create_sql_table(ENTAP_VERSION)) {
                    // YES, create entry into table

                    version_str = get_required_version_str();
                    query = sqlite3_mprintf(
                            "INSERT INTO %Q (%Q) VALUES (%Q);",
                            SQL_TABLE_VERSION_TITLE.c_str(),
                            SQL_TABLE_VERSION_COL_VER.c_str(),
                            version_str.c_str()
                    );

                    FS_dprint("Executing SQL cmd: " + std::string(query));
                    ret = _pDatabaseHelper->execute_cmd(query);

                } else {
                    // NO, return
                    ret = false;
                }
            } else {
                ret = false;
            }
            break;

        default:
            FS_dprint("ERROR: Unhandled set_database_versions type");
            ret = false;
    }
    return ret;
}

