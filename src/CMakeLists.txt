file(GLOB SOURCE_FILES *.cpp expression/*.cpp frame_selection/*.cpp ontology/*.cpp)

include_directories(../libs/pstream)
include_directories(../libs/fast-cpp-csv-parser) # May need to compile with thread control

include_directories(database)
link_directories(database)
add_subdirectory(database)

include_directories(${Boost_INCLUDE_DIRS})
add_executable(EnTAP ${SOURCE_FILES})


target_link_libraries(EnTAP boost pthread EntapDatabase SQLDatabaseHelper sqlite3)
install(TARGETS EnTAP DESTINATION bin)


