cmake_minimum_required(VERSION 3.0)
project(EnTAP VERSION 1.2.1 LANGUAGES CXX C)
include(CheckCXXCompilerFlag)
include(GNUInstallDirs)

option(BUILD_STATIC "BUILD_STATIC" OFF)

if (BUILD_STATIC)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    SET(BUILD_SHARED_LIBRARIES OFF)
    SET(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
    message("Compiler supports C++11!")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()

find_package(PythonInterp)
if (PYTHONINTERP_FOUND)
    message("Python detected! Version: ${PYTHON_VERSION_STRING}")
    message("--   Ensure Python SQLITE and Matplotlib (optional) modules are installed")
else()
    message(WARNING "Python was not detected, EnTAP graphing will not be allowed")
endif()

find_package(Perl)
if (Perl_FOUND)
    message("Perl detected! Version: ${PERL_VERSION_STRING}")
else()
    message(WARNING "Perl was not detected, this is required for GeneMarkS-T!")
endif()

set(SOURCE_FILES
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
        src/database/SQLDatabaseHelper.cpp src/database/SQLDatabaseHelper.h
        src/database/sqlite3.c src/database/sqlite3.h
        src/frame_selection/ModGeneMarkST.cpp src/frame_selection/ModGeneMarkST.h
        src/frame_selection/AbstractFrame.h src/frame_selection/AbstractFrame.cpp
        src/expression/AbstractExpression.h src/expression/AbstractExpression.cpp
        src/expression/ModRSEM.cpp src/expression/ModRSEM.h
        src/ontology/AbstractOntology.h src/ontology/AbstractOntology.cpp
        src/ontology/ModEggnog.cpp src/ontology/ModEggnog.h
        src/ontology/ModInterpro.cpp src/ontology/ModInterpro.h
        src/QueryData.cpp src/QueryData.h src/common.h src/config.h
        src/FileSystem.cpp src/FileSystem.h src/version.h
        src/database/EntapDatabase.cpp src/database/EntapDatabase.h
        src/TerminalCommands.cpp src/TerminalCommands.h
        src/ontology/ModEggnogDMND.h src/ontology/ModEggnogDMND.cpp
        src/database/EggnogDatabase.cpp src/database/EggnogDatabase.h
        src/EntapModule.cpp src/EntapModule.h
        src/similarity_search/AbstractSimilaritySearch.cpp src/similarity_search/AbstractSimilaritySearch.h
        src/similarity_search/ModDiamond.cpp src/similarity_search/ModDiamond.h
        src/QueryAlignment.cpp src/QueryAlignment.h
        src/frame_selection/ModTransdecoder.cpp src/frame_selection/ModTransdecoder.h
        src/database/BuscoDatabase.cpp src/database/BuscoDatabase.h
        src/ontology/ModBUSCO.cpp src/ontology/ModBUSCO.h
        src/tests/UnitTests.cpp src/tests/UnitTests.h
        src/database/NCBIEntrez.cpp src/database/NCBIEntrez.h
        src/user_inputs/GFF_File.cpp src/user_inputs/GFF_File.h src/HorizontalGeneTransfer.cpp
        src/HorizontalGeneTransfer.h src/horizontal_gene_transfer/AbstractHorizontalGeneTransfer.cpp
        src/horizontal_gene_transfer/AbstractHorizontalGeneTransfer.h src/horizontal_gene_transfer/ModHorizontalGeneTransferDiamond.cpp
        src/horizontal_gene_transfer/ModHorizontalGeneTransferDiamond.h
        src/database/NCBIDatabase.cpp
        src/database/NCBIDatabase.h)

# Include libraries
include_directories(libs/pstream)
include_directories(libs/fast-cpp-csv-parser) # May need to compile with thread control
include_directories(libs/cereal-1.2.2/include)
include_directories(libs/tclap-1.2.2/include)
include_directories(libs/boost-1.68.0/include)

set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory(libs)

add_executable(EnTAP ${SOURCE_FILES})

target_link_libraries(EnTAP dl pthread nlohmann_json::nlohmann_json)
install(TARGETS EnTAP DESTINATION bin)