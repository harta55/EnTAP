cmake_minimum_required(VERSION 2.8)
project(enTAP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(SOURCE_FILES libs
        src/main.cpp
        src/ExceptionHandler.cpp src/ExceptionHandler.h
        src/EntapGlobals.h src/EntapGlobals.cpp
        src/EntapConfig.cpp src/EntapConfig.h
        src/EntapExecute.cpp src/EntapExecute.h
        src/QuerySequence.cpp src/QuerySequence.h
        src/FrameSelection.cpp src/FrameSelection.h
        src/ExpressionAnalysis.cpp src/ExpressionAnalysis.h
        src/SimilaritySearch.cpp src/SimilaritySearch.h
        src/Ontology.cpp src/Ontology.h
        src/GraphingManager.cpp src/GraphingManager.h
        src/UserInput.cpp src/UserInput.h
        src/EggnogLevels.h
        src/DatabaseHelper.cpp src/DatabaseHelper.h
        src/sqlite3.c src/sqlite3.h)

add_executable(enTAP ${SOURCE_FILES})

set(Boost_INCLUDE_DIR /home/harta/CLionProjects/enTAP/libs/boost_1_63_0)
set(Boost_LIBRARY_DIR /home/harta/CLionProjects/enTAP/libs/boost_1_63_0/stage/lib)

include_directories(libs/pstream)
include_directories(libs/fast-cpp-csv-parser) # May need to compile with thread control
find_package(Boost 1.63.0 COMPONENTS program_options system filesystem iostreams serialization
        regex REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
        set(BOOST_IOSTREAMS_SOURCE ZLIB)
    include_directories(${ZLIB_INCLUDE_DIRS})
    target_link_libraries(enTAP ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} pthread dl)
endif()