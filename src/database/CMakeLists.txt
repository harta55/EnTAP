set(EntapDatabase_files EntapDatabase.cpp  EntapDatabase.h)
set(SQLDatabaseHelper_files SQLDatabaseHelper.cpp  SQLDatabaseHelper.h)
set(sqlite3_files sqlite3.c  sqlite3.h)

add_library(EntapDatabase ${EntapDatabase_files})
add_library(SQLDatabaseHelper ${SQLDatabaseHelper_files})
add_library(sqlite3 ${sqlite3_files})

target_include_directories(EntapDatabase PUBLIC .)
target_include_directories(SQLDatabaseHelper PUBLIC .)
target_include_directories(sqlite3 PUBLIC .)

target_link_libraries(EntapDatabase boost SQLDatabaseHelper SQLDatabaseHelper)
target_link_libraries(SQLDatabaseHelper sqlite3)
target_link_libraries(sqlite3 dl)

export(TARGETS EntapDatabase SQLDatabaseHelper sqlite3 FILE MyLibraryConfig.cmake)
