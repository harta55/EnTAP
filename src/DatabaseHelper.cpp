//
// Created by harta on 5/25/17.
//

#include "DatabaseHelper.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"

bool DatabaseHelper::open(std::string file) {
    return sqlite3_open(file.c_str(),&_database) == SQLITE_OK;
}

void DatabaseHelper::close() {
    sqlite3_close(_database);
}

std::vector<std::vector<std::string>> DatabaseHelper::query(char *query) {
    sqlite3_stmt *stmt;
    query_struct output;
//    sqlite3_exec(_database,"PRAGMA synchronous = OFF", NULL, NULL, NULL);
//    sqlite3_exec(_database,"PRAGMA count_changes = false", NULL, NULL, NULL);
//    sqlite3_exec(_database,"PRAGMA journal_mode = OFF", NULL, NULL, NULL);
    if (sqlite3_prepare_v2(_database,query,-1,&stmt,0) == SQLITE_OK) {
        int col_num = sqlite3_column_count(stmt);
        int stat = 0;
        while (true) {
            stat = sqlite3_step(stmt);
            if (stat == SQLITE_ROW) {
                std::vector<std::string> vals;
                for (int i = 0; i < col_num; i++) {
                    vals.push_back(std::string((char*)sqlite3_column_text(stmt,i)));
                }
                output.push_back(vals);
            } else {
                break;
            }
        }
        sqlite3_finalize(stmt);
    } else {
        throw ExceptionHandler("Error querying database",ENTAP_ERR::E_INIT_GO_SETUP);
    }
    return output;
}

DatabaseHelper::DatabaseHelper() {
    _database = NULL;
}

DatabaseHelper::~DatabaseHelper() {
    close();
}

