cmake_minimum_required(VERSION 3.6)
project(enTAP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp ExceptionHandler.cpp ExceptionHandler.h ErrorFlags.h EntapInit.cpp EntapInit.h)
add_executable(enTAP ${SOURCE_FILES})

set(Boost_INCLUDE_DIR /home/harta55/ClionProjects/Libraries/boost_1_62_0)
set(Boost_LIBRARY_DIR /home/harta55/ClionProjects/Libraries/boost_1_62_0/stage/lib)
set(CURL_INCLUDE_DIR /home/harta55/CLionProjects/enTAP/libs/curl-7.52.1/include/curl)
set(CURL_LIBRARY /home/harta55/CLionProjects/enTAP/libs/curl-7.52.1/lib)

#set(Boost_INCLUDE_DIR /home/harta/Documents/boost_1_63_0)
#set(Boost_LIBRARY_DIR /home/harta/Documents/boost_1_63_0/stage/lib)
#set(CURL_INCLUDE_DIR /home/harta/CLionProjects/enTAP/libs/curl-7.52.1/include/curl)
#set(CURL_LIBRARY /home/harta/CLionProjects/enTAP/libs/curl-7.52.1/lib)

include_directories(libs/pstream)
find_package(Boost 1.62.0 COMPONENTS program_options REQUIRED)
find_library(CURL_LIBRARY curl curllib libcurl_imp curllib_static)
#find_package(CURL REQUIRED)


if (CURL_FOUND)
    message("CURL HAS BEEN FOUDN!")
    add_definitions(-DCURL_STATICLIB)
endif()



if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    include_directories(${CURL_INCLUDE_DIRS})
    target_link_libraries(enTAP ${Boost_LIBRARIES})
    target_link_libraries(enTAP ${CURL_LIBRARIES})
endif()