cmake_minimum_required(VERSION 2.8)
project(EnTAP)
include(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    message("Compiler supports C++11!")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()

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
        src/ontology/EggnogLevels.h
        src/database/SQLDatabaseHelper.cpp src/database/SQLDatabaseHelper.h
        src/database/sqlite3.c src/database/sqlite3.h
        src/frame_selection/ModGeneMarkST.cpp src/frame_selection/ModGeneMarkST.h
        src/frame_selection/AbstractFrame.h
        src/expression/AbstractExpression.h
        src/expression/ModRSEM.cpp src/expression/ModRSEM.h
        src/ontology/AbstractOntology.h
        src/ontology/ModEggnog.cpp src/ontology/ModEggnog.h
        src/ontology/AbstractOntology.cpp
        src/ontology/ModInterpro.cpp src/ontology/ModInterpro.h
        src/QueryData.cpp src/QueryData.h src/common.h src/config.h
        src/FileSystem.cpp src/FileSystem.h src/version.h
        src/database/EntapDatabase.cpp src/database/EntapDatabase.h
        src/TerminalCommands.cpp src/TerminalCommands.h)

add_executable(EnTAP ${SOURCE_FILES})

# Change these paths to the specific boost directory if you're having trouble linking
set(Boost_INCLUDE_DIR /path/to/boost_1_50_0)
set(Boost_LIBRARY_DIR /path/to/boost_1_50_0/stage/lib)

include_directories(libs/pstream)
include_directories(libs/fast-cpp-csv-parser) # May need to compile with thread control
find_package(Boost 1.50.0 COMPONENTS program_options system filesystem iostreams serialization
        regex date_time REQUIRED)

find_package(PythonInterp)
if (PYTHONINTERP_FOUND)
    message("Python detected! Version: ${PYTHON_VERSION_STRING}")
    message("--   Ensure Python SQLITE and Matplotlib (optional) modules are installed")
else()
    message(WARNING "Python was not detected, EnTAP graphing and EggNOG functionality will not be allowed")
endif()

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(EnTAP ${Boost_LIBRARIES} pthread dl)
else()
    message(FATAL_ERROR "Boost Libraries are required to run EnTAP.")
endif()

target_link_libraries(EnTAP ${Boost_LIBRARIES})
install(TARGETS EnTAP DESTINATION bin)
