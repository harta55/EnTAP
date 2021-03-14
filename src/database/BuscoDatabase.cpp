/******************************************************************
 *
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
 *******************************************************************/

#include "BuscoDatabase.h"
#include "../FileSystem.h"

BuscoDatabase::BuscoDatabase(FileSystem *fileSystem) {
    mpFileSystem = fileSystem;
    mErrCode = ERR_DATA_OK;
    mTempDirectory = mpFileSystem->get_temp_outdir();
}

BuscoDatabase::~BuscoDatabase() {
    FS_dprint("Killing - BuscoDatabase");
}

BuscoDatabase::BUSCO_DB_ERR BuscoDatabase::download_database(std::string &user_input, std::string &output_path) {
    std::string busco_url;
    std::string temp_busco_db;
    mErrCode = ERR_DATA_OK;

    // Get valid database URL
    if (!valid_database(user_input, busco_url)) {
        // ERROR invalid database found
        set_err_msg("ERROR invalid BUSCO database input (" + user_input + ") and URL could not be found",
                    ERR_DATA_UNKNOWN_DATABASE);
        return ERR_DATA_UNKNOWN_DATABASE;   // EXIT - ERROR unknown database
    }

    // Have a valid database URL, download
    temp_busco_db = PATHS(mTempDirectory, get_database_shortname(busco_url));
    if (!mpFileSystem->download_ftp_file(busco_url, temp_busco_db)) {
        // ERROR could not download database
        set_err_msg("ERROR unable to download BUSCO database from: " + busco_url +
                    mpFileSystem->get_error(), ERR_DATA_DOWNLOAD);
        return ERR_DATA_DOWNLOAD;       // EXIT - ERROR unable to download database
    }

    // Decompress/move the file
    if (!mpFileSystem->decompress_file(temp_busco_db, output_path, FileSystem::ENT_FILE_TAR_GZ)) {
        // ERROR could not decompress/move the database
        set_err_msg("ERROR unable to decompress the BUSCO database at: " + temp_busco_db +
                    mpFileSystem->get_error(), ERR_DATA_DECOMPRESS);
        return ERR_DATA_DECOMPRESS;     // EXIT - ERROR unable to decompress/move
    }

    // remove temporary file, ignore errors
    mpFileSystem->delete_dir(temp_busco_db);

    return ERR_DATA_OK;
}

bool BuscoDatabase::valid_database(std::string &database, std::string &database_url) {
    bool ret = false;           // TRUE if valid URL found, FALSE if otherwise
    std::string ret_str="";     // Returned database_url

        // Is this a database we already have a mapping for?
    auto it = DATASET_TO_URL.find(database);
    if (it != DATASET_TO_URL.end()) {
        // Yes, use the URL we already have mapped
        ret = true;
        ret_str = it->second;
    } else {
        // No, check whether the URL is valid
        if (mpFileSystem->is_url(database)) {
            ret = true;
            ret_str = database;
        } else {
            ; // ERROR invalid database
        }
    }
    database_url = ret_str;
    return ret;
}

std::string BuscoDatabase::print_error_log() {
    return "\nBUSCO Database Error: " + mErrMsg;
}

void BuscoDatabase::set_err_msg(std::string msg, BuscoDatabase::BUSCO_DB_ERR code) {
    FS_dprint(msg);
    mErrCode = code;
    mErrMsg = msg;
}

// Returns database shortname from either url input or general database
std::string BuscoDatabase::get_database_shortname(std::string &database) {
    std::string ret="";
    std::string temp;

    if (valid_database(database, temp)) {
        if (mpFileSystem->is_url(database)) {
            ret = mpFileSystem->get_filename(database, false);
        } else {
            ret = database;
        }
    } else {
        ; // ERROR Invalid database
    }
    return ret;
}
