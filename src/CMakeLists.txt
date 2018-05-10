set(SOURCE_FILES ../libs
        main.cpp
        ExceptionHandler.cpp ExceptionHandler.h
        EntapGlobals.h EntapGlobals.cpp
        EntapConfig.cpp EntapConfig.h
        EntapExecute.cpp EntapExecute.h
        QuerySequence.cpp QuerySequence.h
        FrameSelection.cpp FrameSelection.h
        ExpressionAnalysis.cpp ExpressionAnalysis.h
        SimilaritySearch.cpp SimilaritySearch.h
        Ontology.cpp Ontology.h
        GraphingManager.cpp GraphingManager.h
        UserInput.cpp UserInput.h
        ontology/EggnogLevels.h
        frame_selection/ModGeneMarkST.cpp frame_selection/ModGeneMarkST.h
        frame_selection/AbstractFrame.h
        expression/AbstractExpression.h
        expression/ModRSEM.cpp expression/ModRSEM.h
        ontology/AbstractOntology.h
        ontology/ModEggnog.cpp ontology/ModEggnog.h
        ontology/AbstractOntology.cpp
        ontology/ModInterpro.cpp ontology/ModInterpro.h
        QueryData.cpp QueryData.h common.h config.h
        FileSystem.cpp FileSystem.h version.h
        TerminalCommands.cpp TerminalCommands.h)

include_directories(../libs/pstream)
include_directories(../libs/fast-cpp-csv-parser) # May need to compile with thread control

include_directories(database)
link_directories(database)
add_subdirectory(database)

include_directories(${Boost_INCLUDE_DIRS})
add_executable(EnTAP ${SOURCE_FILES})


target_link_libraries(EnTAP boost pthread EntapDatabase SQLDatabaseHelper sqlite3)
install(TARGETS EnTAP DESTINATION bin)


