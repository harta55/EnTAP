/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
 *
 * This file is part of EnTAP.
 *
 * EnTAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EnTAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "FileSystem.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
#include <iomanip>
#include <iostream>
#include <sys/stat.h>
#include <chrono>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/time_clock.hpp>
#include "config.h"
#include <ctime>
#include <cstring>
#include <boost/date_time/posix_time/posix_time.hpp>

const std::string FileSystem::EXT_TXT = ".txt";
const std::string FileSystem::EXT_ERR = ".err";
const std::string FileSystem::EXT_OUT = ".out";
const std::string FileSystem::EXT_BAM = ".bam";
const std::string FileSystem::EXT_FAA = ".faa";
const std::string FileSystem::EXT_FNN = ".fnn";
const std::string FileSystem::EXT_XML = ".xml";
const std::string FileSystem::EXT_DMND= ".dmnd";


// Removed for older compilers, may bring back
void FileSystem::open_out(std::string &path, std::ofstream &ofstream) {
    ofstream = std::ofstream(path,std::ios::out | std::ios::app);
    if (!ofstream.is_open()) {
        throw ExceptionHandler("Error opening file: " + path,
             ERR_ENTAP_FILE_IO);
    }
}


/**
 * ======================================================================
 * Function void FS_close_file(std::ofstream &ofstream)
 *
 * Description          - Close a stream
 *
 * Notes                - None
 *
 * @param ofstream      - File stream
 *
 * @return              - None
 *
 * =====================================================================
 */
void FileSystem::close_file(std::ofstream &ofstream) {
    try {
        ofstream.close();
    } catch (const std::exception &exception) {
        throw ExceptionHandler(exception.what(), ERR_ENTAP_FILE_IO);
    }
}


/**
 * ======================================================================
 * Function void print_debug(std::string    msg)
 *
 * Description          - Handles printing to EnTAP debug file
 *                      - Adds timestamp to each entry
 *
 * Notes                - None
 *
 * @param msg           - Message to be sent to debug file
 * @return              - None
 *
 * =====================================================================
 */
void FS_dprint(const std::string &msg) {

#if DEBUG
    std::chrono::time_point<std::chrono::system_clock> current;
    std::time_t time;

    current = std::chrono::system_clock::now();
    time = std::chrono::system_clock::to_time_t(current);
    std::string out_time(std::ctime(&time));
    std::ofstream debug_file(DEBUG_FILE_PATH, std::ios::out | std::ios::app);

    debug_file << out_time.substr(0,out_time.length()-1) << ": " + msg << std::endl;
    debug_file.close();
#endif
}


/**
 * ======================================================================
 * Function void print_statistics(std::string    &msg)
 *
 * Description          - Handles printing to EnTAP statistics/log file
 *
 * Notes                - None
 *
 * @param msg           - Message to be sent to log file
 * @return              - None
 *
 * =====================================================================
 */
void FileSystem::print_stats(std::string &msg) {
    std::ofstream log_file(LOG_FILE_PATH, std::ios::out | std::ios::app);
    log_file << msg << std::endl;
    close_file(log_file);
}


/**
 * ======================================================================
 * Function bool file_exists(std::string path)
 *
 * Description          - Checks whether a file exists in the OS
 *                      - Currently multiplatform boost implementation
 *
 * Notes                - None
 *
 * @param path          - Path of file to verify
 *
 * @return              - true/false of file existing
 *
 * =====================================================================
 */
bool FileSystem::file_exists(std::string path) {
#ifdef USE_BOOST
    return boost::filesystem::exists(path);
#else
    struct stat buff;
    return (stat(path.c_str(), &buff) == 0);
#endif
}


/**
 * ======================================================================
 * Function bool FS_file_is_open(std::ofstream &ofstream)
 *
 * Description          - Checks if an ofstream is opened
 *
 * Notes                - Unused b/c older compilers
 *
 * @param ofstream      - Stream
 *
 * @return              - True/false if open
 * ======================================================================
 */
bool FileSystem::file_is_open(std::ofstream &ofstream) {
    return ofstream.is_open();
}


/**
 * ======================================================================
 * Function bool FS_file_test_open(std::string &path)
 *
 * Description          - Checks whether a file can be opened successfully
 *
 * Notes                - None
 *
 * @param path          - Path to file
 *
 * @return              - True/false if successful
 * ======================================================================
 */
bool FileSystem::file_test_open(std::string &path) {
    bool is_open;
    std::ifstream ifstream(path);
    is_open = ifstream.is_open();
    ifstream.close();
    return is_open;
}


/**
 * ======================================================================
 * Function bool FS_delete_file(std::string path)
 *
 * Description          - Delete target file
 *
 * Notes                - None
 *
 * @param path          - Path to file
 *
 * @return              - True/false if successful
 * ======================================================================
 */
bool FileSystem::delete_file(std::string path) {
    if (!file_exists(path)) return false;
#ifdef USE_BOOST
    FS_dprint("Deleting file: " + path);
    return boostFS::remove(path);
#else
    return remove(path) == 0;
#endif
}


/**
 * ======================================================================
 * Function bool FS_directory_iterate(bool del, std::string &path)
 *
 * Description          - Boost recursive iteration through directory
 *
 * Notes                - None
 *
 * @param path          - Path to directory
 * @param del           - True/false to delete files that are empty (unused)
 *
 * @return              - True/false if successful
 * ======================================================================
 */
bool FileSystem::directory_iterate(bool del, std::string &path) {
    FS_dprint("Iterating through directory: " + path);
    if (!file_exists(path)) return false;
    try {
        for (boostFS::recursive_directory_iterator it(path), end; it != end; ++it) {
            if (!boostFS::is_directory(it->path())) {
                // Is file
                if (file_empty(it->path().string()) && del) {
                    delete_file(it->path().string());
                    FS_dprint("Deleted: " + it->path().string());
                }
            }
        }
    } catch (...) {
        return false;
    }
    FS_dprint("Success!");
    return true;
}


/**
 * ======================================================================
 * Function bool FS_file_empty(std::string path)
 *
 * Description          - Check if a file is empty
 *                      - Checks no line or all empty lines
 *
 * Notes                - None
 *
 * @param path          - Path to file
 *
 * @return              - True/false if file empty
 * ======================================================================
 */
bool FileSystem::file_empty(std::string path) {
    std::ifstream file(path);
    bool empty;
    std::string line;
    empty = file.peek() == std::ifstream::traits_type::eof();
    if (!empty) {
        empty = true;
        while (getline(file,line)) {
            if (line.empty()) continue;
            empty = false;
            break;
        }
    }
    file.close();
    return empty;
}


/**
 * ======================================================================
 * Function bool FS_check_fasta(std::string& path)
 *
 * Description          - Minor check on fasta file for format
 *
 * Notes                - None
 *
 * @param path          - Path to fasta file
 *
 * @return              - True/false if fasta file valid
 * ======================================================================
 */
bool FileSystem::check_fasta(std::string& path) {
    std::string line;
    bool valid = false;
    try {
        std::ifstream file(path);
        while (getline(file, line)) {
            if (line.at(0) == '>') {
				valid = true;
				break;
			}
        }
        file.close();
        return valid;
    } catch (...) {
        return false;
    }
}


/**
 * ======================================================================
 * Function bool FS_create_dir(std::string& path)
 *
 * Description          - Create directory
 *
 * Notes                - None
 *
 * @param path          - Path to directory
 *
 * @return              - True/false if created successfully
 * ======================================================================
 */
bool FileSystem::create_dir(std::string& path) {
#ifdef USE_BOOST
    return boostFS::create_directories(path);
#endif
}


/**
 * ======================================================================
 * Function bool FS_delete_dir(std::string& path)
 *
 * Description          - Delete directory
 *
 * Notes                - None
 *
 * @param path          - Path to directory
 *
 * @return              - None
 * ======================================================================
 */
void FileSystem::delete_dir(std::string& path) {
#ifdef USE_BOOST
    boostFS::remove_all(path);
#endif
}


/**
 * ======================================================================
 * Function FS_file_no_lines(std::string path)
 *
 * Description          - Check if specific file is empty (has no lines)
 *
 * Notes                - Used for certain RSEM execution that does not
 *                        relay error code of non-zero on failure
 *
 * @param path          - Path to file
 *
 * @return              - True if file is empty
 *
 * =====================================================================
 */
bool FileSystem::file_no_lines(std::string path) {
    std::ifstream ifstream(path);
    return ifstream.peek() == std::ifstream::traits_type::eof();
}


/**
 * ======================================================================
 * Function std::vector<std::string> FS_list_to_vect(char it, std::string &list)
 *
 * Description          - Not FS but just threw here...
 *                      - Turns string of list (red,blue,green) to a vector of these
 *
 * Notes                - None
 *
 * @param it            - Iterating character
 * @param list          - String of list (ie: red,blue,green)
 *
 * @return              - True if file is empty
 *
 * =====================================================================
 */
std::vector<std::string> FileSystem::list_to_vect(char it, std::string &list) {
    std::string temp;
    std::vector<std::string> output;

    if (list.empty()) return std::vector<std::string>();

    std::istringstream vals(list);
    while (std::getline(vals,temp,it)) {
        output.push_back(temp);
    }
    return output;
}


/**
 * ======================================================================
 * Function std::string FS_get_cur_dir()
 *
 * Description          - Returns the current working directory of execution
 *
 * Notes                - None
 *
 * @return              - Current directory
 *
 * =====================================================================
 */
std::string FileSystem::get_cur_dir() {

#ifdef USE_BOOST
    return boostFS::current_path().string();
#endif

}

FileSystem::~FileSystem() {
    FS_dprint("Killing Object - FileSystem");

}

FileSystem::FileSystem(std::string &root) {
    // This routine will process entire root directory here and generate
    // hierarchy

    _root_path = root;
    create_dir(root);
    init_log();
}

/**
 * ======================================================================
 * Function init_log()
 *
 * Description          - Initializes log and debug files
 *
 * Notes                - None
 *
 * @return              - None
 * ======================================================================
 */
void FileSystem::init_log() {
#ifndef USE_BOOST
    std::chrono::time_point<std::chrono::system_clock> now;
    std::time_t                                        now_time;
    std::tm                                            now_tm;
#endif
    std::stringstream                                  ss;
    std::string                                        log_file_name;
    std::string                                        debug_file_name;
    std::string                                        time_date;

#ifdef USE_BOOST
    boost::posix_time::ptime local = boost::posix_time::second_clock::local_time();
    ss <<
       "_" <<
       local.date().year()             << "." <<
       local.date().month().as_number()<< "." <<
       local.date().day()              << "-" <<
       local.time_of_day().hours()     << "h"   <<
       local.time_of_day().minutes()   << "m" <<
       local.time_of_day().seconds()   << "s";
#else
    now      = std::chrono::system_clock::now();
    now_time = std::chrono::system_clock::to_time_t(now);
    now_tm   = *std::localtime(&now_time);
    ss << std::put_time(&now_tm, "_%Y.%m.%d-%Hh.%Mm.%Ss");
#endif
    time_date       = ss.str();
    log_file_name   = LOG_FILENAME   + time_date + LOG_EXTENSION;
    debug_file_name = DEBUG_FILENAME + time_date + LOG_EXTENSION;
    DEBUG_FILE_PATH = PATHS(_root_path, debug_file_name);
    LOG_FILE_PATH   = PATHS(_root_path, log_file_name);
    delete_file(DEBUG_FILE_PATH);
    delete_file(LOG_FILE_PATH);
    FS_dprint("Start - EnTAP");
}

const std::string &FileSystem::get_root_path() const {
    return _root_path;
}

std::string FileSystem::get_file_extension(const std::string &path) {
#ifdef USE_BOOST
    boostFS::path bpath(path);
    return bpath.extension().string();
#endif
}

bool FileSystem::copy_file(std::string inpath, std::string outpath, bool overwrite) {
#ifdef USE_BOOST
    try {
        if (overwrite) {
            boostFS::copy_file(inpath,outpath,boostFS::copy_option::overwrite_if_exists);
        } else {
            boostFS::copy_file(inpath,outpath);
        }
    } catch (...) {
        return false;
    }
    return true;
#endif
}
