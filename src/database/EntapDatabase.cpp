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

#include "EntapDatabase.h"

EntapDatabase::EntapDatabase(FileSystem* filesystem) {
    // Initialize
    _pFilesystem     = filesystem;
    _temp_directory  = filesystem->get_temp_outdir();    // created previously
}

EntapDatabase::DATABASE_ERR EntapDatabase::set_database(DATABASE_TYPE type) {

}

EntapDatabase::DATABASE_ERR EntapDatabase::download_database(EntapDatabase::DATABASE_TYPE type, std::string &path) {
    switch (type) {
        case ENTAP_SQL:
            download_entap_sql(path);
            break;
        default:
            return DATA_OK;
    }
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_database(EntapDatabase::DATABASE_TYPE type, std::string &path) {
    switch (type) {
        case ENTAP_SQL:
            download_entap_sql(path);
        default:
            return DATA_OK;
    }
}

EntapDatabase::DATABASE_ERR EntapDatabase::download_entap_sql(std::string &outpath) {
    DATABASE_ERR err_code;

    FS_dprint("Creating EnTAP SQL database...");

    if (_pDatabaseHelper != nullptr) {
        // SQL should not already have been created
        return DATA_SQL_DUPLICATE;
    } else if (_pFilesystem->file_exists(outpath)) {
        // Out file should NOT exist yet
        return DATA_FILE_EXISTS;
    }

    // Create sql database (will create if doesn't exist)
    FS_dprint("Creating SQL database...");
    _pDatabaseHelper = new SQLDatabaseHelper();
    if (!_pDatabaseHelper->create(outpath)) {
        // Return error if can't create
        return DATA_SQL_CREATE;
    }
    FS_dprint("Success!");

    // Generate tax entries, don't need a path - using SQL member
    err_code = generate_entap_tax(ENTAP_SQL, "");
    if (err_code != DATA_OK) return err_code;

    // Generate go entries







    return DATA_OK;
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_sql(std::string &) {
    return DATA_OK;
}

EntapDatabase::~EntapDatabase() {
    if (_pDatabaseHelper != nullptr) {
        _pDatabaseHelper->close();
        delete(_pDatabaseHelper);
    }
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_tax(EntapDatabase::DATABASE_TYPE, std::string) {
    std::string temp_outpath;   // Path

    FS_dprint("Generating EnTAP Tax database entries");

    return DATA_OK;
}

EntapDatabase::DATABASE_ERR EntapDatabase::generate_entap_go(EntapDatabase::DATABASE_TYPE, std::string) {
    return DATA_OK;
}


