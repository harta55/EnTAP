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


EntapDatabase::EntapDatabase(FileSystem* filesystem) {
    // Initialize
    _pFilesystem     = filesystem;
    _temp_directory  = filesystem->get_temp_outdir();    // created previously
    _pSerializedDatabase = nullptr;
    _pDatabaseHelper     = nullptr;
}

EntapDatabase::DATABASE_ERR EntapDatabase::set_database(DATABASE_TYPE type) {
    return ERR_DATA_OK;
}

EntapDatabase::DATABASE_ERR EntapDatabase::download_database(EntapDatabase::DATABASE_TYPE type, std::string &path) {
    switch (type) {
        case ENTAP_SQL:
            download_entap_sql(path);
            break;
        default:
            return ERR_DATA_OK;
    }
    return ERR_DATA_OK;
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_database(EntapDatabase::DATABASE_TYPE type, std::string &path) {
    DATABASE_ERR err;

    switch (type) {
        case ENTAP_SQL:
            err = generate_entap_sql(path);
            break;
        case ENTAP_SERIALIZED:
            err = generate_entap_serial(path);
            break;
        default:
            return ERR_DATA_OK;
    }
    if (err != ERR_DATA_OK) _pFilesystem->delete_file(path);
    return err;
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_sql(std::string &outpath) {
    DATABASE_ERR err_code;

    FS_dprint("Creating EnTAP SQL database...");

    if (_pDatabaseHelper != nullptr) {
        // SQL should not already have been created
        return ERR_DATA_SQL_DUPLICATE;
    } else if (_pFilesystem->file_exists(outpath)) {
        // Out file should NOT exist yet
        return ERR_DATA_FILE_EXISTS;
    }

    // Create sql database (will create if doesn't exist)
    FS_dprint("Creating SQL database...");
    _pDatabaseHelper = new SQLDatabaseHelper();
    if (!_pDatabaseHelper->create(outpath)) {
        // Return error if can't create
        return ERR_DATA_SQL_CREATE_DATABASE;
    }
    FS_dprint("Success!");

    // Generate tax entries, don't need a path - using SQL member
    err_code = generate_entap_tax(ENTAP_SQL, "");
    if (err_code != ERR_DATA_OK) {
        return err_code;
    }

    // Generate go entries
    err_code = generate_entap_go(ENTAP_SQL, "");
    if (err_code != ERR_DATA_OK) {
        return err_code;
    }

    return ERR_DATA_OK;
}

EntapDatabase::DATABASE_ERR EntapDatabase::download_entap_sql(std::string &) {
    return ERR_DATA_OK;
}

EntapDatabase::~EntapDatabase() {
    if (_pDatabaseHelper != nullptr) {
        _pDatabaseHelper->close();
        delete(_pDatabaseHelper);
    }
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_tax(EntapDatabase::DATABASE_TYPE type,
                                                              std::string outpath) {
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
            FS_dprint("Error generating SQL taxonomy table");
            return ERR_DATA_SQL_CREATE_TABLE;
        }
    }

    temp_outpath = PATHS(_temp_directory, NCBI_TAX_DUMP_FILENAME);

    // download files from NCBI tax
    if (!_pFilesystem->download_ftp_file(FTP_NCBI_TAX_DUMP_TARGZ, temp_outpath)) {
        return ERR_DATA_TAX_DOWNLOAD;
    }
    // decompress TAR.GZ file
    if (!_pFilesystem->decompress_file(temp_outpath, _temp_directory, FileSystem::FILE_TAR_GZ)) {
        return ERR_DATA_FILE_DECOMPRESS;
    }
    // remove tar file
    _pFilesystem->delete_file(temp_outpath);
    // set paths to uncompressed files
    ncbi_names_path = PATHS(_temp_directory, NCBI_TAX_DUMP_FTP_NAMES);
    ncbi_nodes_path = PATHS(_temp_directory, NCBI_TAX_DUMP_FTP_NODES);

    FS_dprint("Files downloaded and compressed, parsing...");


    FS_dprint("Parsing NCBI Names file at: " + ncbi_names_path);
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
            current_entries++;
            taxEntry = {};
            taxEntry.lineage = entap_tax_get_lineage(pair.second, taxonomy_nodes);
            taxEntry.tax_id  = pair.second.ncbi_id;
            taxEntry.tax_name= name;

            // Add to SQL database or other...
            if (type == ENTAP_SQL) {
                if (!sql_add_tax_entry(taxEntry)) {
                    // unable to add entry
                    FS_dprint("Unable to add tax entry: " + name);
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
    FS_dprint("Success! NCBI data complete");
    return ERR_DATA_OK;
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_go(EntapDatabase::DATABASE_TYPE type,
                                                             std::string outpath) {
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
        return ERR_DATA_GO_DOWNLOAD;
    }

    // decompress database file
    if (!_pFilesystem->decompress_file(go_database_targz, _temp_directory, FileSystem::FILE_TAR_GZ)) {
        // failed to decompress
        return ERR_DATA_GO_DECOMPRESS;
    }
    _pFilesystem->delete_file(go_database_targz);

    // Files are packaged within a directory. Set paths
    go_database_dir = PATHS(_temp_directory, GO_TERMDB_DIR);
    go_term_path    = PATHS(go_database_dir, GO_TERM_FILE);
    go_graph_path   = PATHS(go_database_dir, GO_GRAPH_FILE);

    if (!_pFilesystem->file_exists(go_term_path) || !_pFilesystem->file_exists(go_graph_path)) {
        FS_dprint("Necessary Gene Ontology files do not exist at:\n" + go_term_path +
                "\n" + go_graph_path);
        return ERR_DATA_GO_DOWNLOAD;
    }

    // If we are creating SQL database, add GO table
    if (type == ENTAP_SQL) {
        if (!create_sql_table(ENTAP_GENE_ONTOLOGY)) {
            // error creating table
            FS_dprint("Unable to create Gene Ontology SQL table");
            return ERR_DATA_SQL_CREATE_TABLE;
        }
    }

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
                FS_dprint("Unable to add GO entry: " + goEntry.go_id);
                return ERR_DATA_GO_ENTRY;
            }
        } else {
            _pSerializedDatabase->gene_ontology_data[go] = goEntry;
        }
    }
    FS_dprint("Success! Gene Ontology data complete");

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

bool EntapDatabase::create_sql_table(DATABASE_TYPE type) {
    std::string table_title;
    char *sql_cmd;

    if (_pDatabaseHelper == nullptr) return false;
    FS_dprint("Creating table...");

    switch (type) {
        case ENTAP_TAXONOMY:
            table_title = SQL_TABLE_NCBI_TAX_TITLE;
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
            return _pDatabaseHelper->execute_cmd(sql_cmd);

        case ENTAP_GENE_ONTOLOGY:
            table_title = SQL_TABLE_GO_TITLE;
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
            return _pDatabaseHelper->execute_cmd(sql_cmd);
        default:
            return false;
    }
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

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_serial(std::string &) {
    return ERR_DATA_OK;
}

EntapDatabase::TaxonomyNode::TaxonomyNode(std::string id) {
    parent_id = "";
    sci_name = "";
    ncbi_id = id;
}