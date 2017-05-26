cmake_minimum_required(VERSION 2.8)
project(enTAP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES libs main.cpp ExceptionHandler.cpp ExceptionHandler.h EntapConsts.h
        EntapInit.cpp EntapInit.h EntapExecute.cpp EntapExecute.h QuerySequence.cpp
        QuerySequence.h FrameSelection.cpp FrameSelection.h ExpressionAnalysis.cpp
        ExpressionAnalysis.h SimilaritySearch.cpp SimilaritySearch.h Ontology.cpp
        Ontology.h sqlite3.h sqlite3.c DatabaseHelper.cpp DatabaseHelper.h)

add_executable(enTAP ${SOURCE_FILES})

#set(Boost_INCLUDE_DIR /home/harta55/ClionProjects/Libraries/boost_1_62_0)
#set(Boost_LIBRARY_DIR /home/harta55/ClionProjects/Libraries/boost_1_62_0/stage/lib)

set(Boost_INCLUDE_DIR /home/harta/Documents/boost_1_63_0)
set(Boost_LIBRARY_DIR /home/harta/Documents/boost_1_63_0/stage/lib)

include_directories(libs/sqlite)
include_directories(libs/pstream)
include_directories(libs/fast-cpp-csv-parser) # May need to compile with thread control
find_package(Boost 1.62.0 COMPONENTS program_options system filesystem iostreams serialization
        regex REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
        set(BOOST_IOSTREAMS_SOURCE ZLIB)
    include_directories(${ZLIB_INCLUDE_DIRS})
    target_link_libraries(enTAP ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} pthread dl)
endif()