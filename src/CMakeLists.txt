file(GLOB SOURCE_FILES *.cpp *.h */*.h */*.cpp */*.c)

include_directories(../libs/pstream)
include_directories(../libs/fast-cpp-csv-parser) # May need to compile with thread control

include_directories(${Boost_INCLUDE_DIRS})
add_executable(EnTAP ${SOURCE_FILES})

target_link_libraries(EnTAP dl pthread Boost::program_options Boost::system Boost::filesystem Boost::iostreams Boost::serialization Boost::regex Boost::date_time)
install(TARGETS EnTAP DESTINATION bin)


