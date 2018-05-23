cmake_minimum_required(VERSION 3.0)
project(EnTAP VERSION 0.8.2 LANGUAGES CXX C)
include(CheckCXXCompilerFlag)
include(GNUInstallDirs)

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
    message("Compiler supports C++11!")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()

option(Boost_USE_STATIC_LIBS "Link with static Boost libraries" OFF)
if(Boost_USE_STATIC_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
else()
    add_definitions( -DBOOST_ALL_DYN_LINK )
endif()

set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.50.0 COMPONENTS program_options system filesystem iostreams serialization regex date_time REQUIRED)
if(Boost_FOUND)
    add_library(boost INTERFACE IMPORTED)
    set_property(TARGET boost PROPERTY
       INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS})
    set_property(TARGET boost PROPERTY
       INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
else()
    message(FATAL_ERROR "Boost Libraries are required to run EnTAP.")
endif()

find_package(PythonInterp)
if (PYTHONINTERP_FOUND)
    message("Python detected! Version: ${PYTHON_VERSION_STRING}")
    message("--   Ensure Python SQLITE and Matplotlib (optional) modules are installed")
else()
    message(WARNING "Python was not detected, EnTAP graphing and EggNOG functionality will not be allowed")
endif()

add_subdirectory(src)
