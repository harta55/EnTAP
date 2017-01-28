cmake_minimum_required(VERSION 3.6)
project(enTAP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp ExceptionHandler.cpp ExceptionHandler.h)
#add_executable(enTAP ${SOURCE_FILES})

#set(Boost_INCLUDE_DIR /home/harta55/ClionProjects/Libraries/boost_1_62_0)
#set(Boost_LIBRARY_DIR /home/harta55/ClionProjects/Libraries/boost_1_62_0/stage/lib)

set(Boost_INCLUDE_DIR /home/harta/Documents/boost_1_63_0)
set(Boost_LIBRARY_DIR /home/harta/Documents/boost_1_63_0/stage/lib)



include_directories(libs/pstream)
find_package(Boost 1.62.0 COMPONENTS program_options REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(enTAP ${SOURCE_FILES})
    target_link_libraries(enTAP ${Boost_LIBRARIES})
endif()