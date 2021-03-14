/*
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
*/


//*********************** Includes *****************************
#include "SQLDatabaseHelper.h"
#include "../ExceptionHandler.h"
#include "../EntapGlobals.h"
//**************************************************************


/**
 * ======================================================================
 * Function bool DatabaseHelper::open(std::string file)
 *
 * Description          - Opens sql database through sqlite3
 *
 * Notes                - This is used as a creation routine as well
 *
 * @param file          - Path to database
 *
 * @return              - True/false if successful
 *
 * =====================================================================
 */
bool SQLDatabaseHelper::open(std::string file) {
    FS_dprint("Opening SQL database at: " + file);
    int err_code;
    err_code = sqlite3_open(file.c_str(),&mpDatabase);
    if (err_code == SQLITE_OK) {
        sqlite3_exec(mpDatabase,"PRAGMA synchronous = OFF", NULL, NULL, NULL);
        sqlite3_exec(mpDatabase,"PRAGMA count_changes = false", NULL, NULL, NULL);
        sqlite3_exec(mpDatabase,"PRAGMA journal_mode = OFF", NULL, NULL, NULL);
        FS_dprint("Success!");
    }
    return err_code == SQLITE_OK;
}


/**
 * ======================================================================
 * Function void DatabaseHelper::close()
 *
 * Description          - Close SQL database
 *
 * Notes                - None
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void SQLDatabaseHelper::close() {
    sqlite3_close(mpDatabase);
}


/**
 * ======================================================================
 * Function bool DatabaseHelper::create(std::string file)
 *
 * Description          - Creates SQL database
 *                      - Calls open (which uses SQLITE_OPEN_CREATE)
 *
 * Notes                - None
 *
 * @param file          - Path to database
 *
 * @return              - True/false if successful
 *
 * =====================================================================
 */
bool SQLDatabaseHelper::create(std::string file) {
    return open(file);
}


/**
 * ======================================================================
 * Function std::vector<std::vector<std::string>> DatabaseHelper::query(char *query)
 *
 * Description          - Queries SQL database and returns vector of
 *                        relevant information pulled from database
 *
 * Notes                - None
 *
 * @param query         - SQL query
 *
 * @return              - Vector of queried information
 *
 * =====================================================================
 */
std::vector<std::vector<std::string>> SQLDatabaseHelper::query(char *query) {
    sqlite3_stmt *stmt;
    query_struct output;
    char* txt;
    if (sqlite3_prepare_v2(mpDatabase,query,-1,&stmt,0) == SQLITE_OK) {
        int col_num = sqlite3_column_count(stmt);
        int stat = 0;
        while (true) {
            stat = sqlite3_step(stmt);
            if (stat == SQLITE_ROW) {
                std::vector<std::string> vals;
                for (int i = 0; i < col_num; i++) {
                    try {
                        txt = (char*)sqlite3_column_text(stmt, i);
                        if (txt != nullptr) {
                            vals.push_back(std::string(txt));
                        } else {
                            vals.push_back("");
                        }
                    } catch (const std::exception &e) {
                        FS_dprint("SQL Error: " + std::string(e.what()));
                    }
                }
                output.push_back(vals);
            } else {
                break;
            }
        }
        sqlite3_finalize(stmt);
    } else {
        throw ExceptionHandler("Error querying database: " + std::string(sqlite3_errmsg(mpDatabase)),
                               ERR_ENTAP_DATABASE_QUERY);
    }
    return output;
}


SQLDatabaseHelper::SQLDatabaseHelper() {
    mpDatabase = NULL;
}


SQLDatabaseHelper::~SQLDatabaseHelper() {
    close();
}

bool SQLDatabaseHelper::execute_cmd(char* cmd) {
    int err;
    char* err_msg = 0;

//    FS_dprint("Executing SQL cmd:\n" + cmd);

    err = sqlite3_exec(mpDatabase, cmd, NULL, 0, &err_msg);

    if( err != SQLITE_OK ){
        std::string err_str = err_msg;
        FS_dprint("SQL Error: " + err_str);
        sqlite3_free(err_msg);
        return false;
    } else {
//        FS_dprint("Success!");
        return true;
    }
}

std::string SQLDatabaseHelper::format_container(std::set<std::string> &in_cont) {
    std::string ret = "(";

    if (in_cont.empty()) return "";
    for (const std::string &val : in_cont) {
        ret += "'" + val + "',";
    }
    ret.erase(ret.find_last_of(','));
    ret += ")";
    return ret;
}

std::string SQLDatabaseHelper::format_string(std::string &str, char delim) {
    std::string ret = "(";

    // Remove newline if exists
    if (str.empty()) return "";
    str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
    std::istringstream iss(str);
    std::string val;
    while(std::getline(iss, val, delim)) {
        ret += "'" + val + "',";
    }
    ret.erase(ret.find_last_of(','));
    ret += ")";
    return ret;
}
