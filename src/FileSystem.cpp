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
#include "config.h"

#if 0 // Removed for older compilers, will bring back later
void FS_open_out(std::string &path, std::ofstream &ofstream) {
    ofstream = std::ofstream(path,std::ios::out | std::ios::app);
    if (!ofstream.is_open()) {
        throw ExceptionHandler("Error opening file: " + path,
            ENTAP_ERR::E_FILE_IO);
    }
}
#endif

void FS_close_file(std::ofstream &ofstream) {
    try {
        ofstream.close();
    } catch (const std::exception &exception) {
        throw ExceptionHandler(exception.what(), ENTAP_ERR::E_FILE_IO);
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
void FS_dprint(std::string msg) {

#if DEBUG
    std::chrono::time_point<std::chrono::system_clock> current;
    std::time_t time;

    current = std::chrono::system_clock::now();
    time = std::chrono::system_clock::to_time_t(current);
    std::string out_time(std::ctime(&time));
    std::ofstream debug_file(DEBUG_FILE_PATH, std::ios::out | std::ios::app);

    debug_file << out_time.substr(0,out_time.length()-1) << ": " + msg << std::endl;
    FS_close_file(debug_file);
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
void FS_print_stats(std::string &msg) {
    std::ofstream log_file(LOG_FILE_PATH, std::ios::out | std::ios::app);
    log_file << msg << std::endl;
    FS_close_file(log_file);
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
bool FS_file_exists(std::string path) {
#ifdef USE_BOOST
    return boost::filesystem::exists(path);
#else
    struct stat buff;
    return (stat(path.c_str(), &buff) == 0);
#endif
}


bool FS_file_is_open(std::ofstream &ofstream) {
    return ofstream.is_open();
}

bool FS_file_test_open(std::string &path) {
    bool is_open;
    std::ifstream ifstream(path);
    is_open = ifstream.is_open();
    ifstream.close();
    return is_open;
}

bool FS_delete_file(std::string path) {
    FS_dprint("Deleting file: " + path);
    if (!FS_file_exists(path)) return false;
    return boostFS::remove(path);
}

bool FS_directory_iterate(bool del, std::string &path) {
    FS_dprint("Iterating through directory: " + path);
    if (!FS_file_exists(path)) return false;
    try {
        for (boostFS::recursive_directory_iterator it(path), end; it != end; ++it) {
            if (!boostFS::is_directory(it->path())) {
                // Is file
                if (FS_file_empty(it->path().string()) && del) {
                    FS_delete_file(it->path().string());
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


bool FS_file_empty(std::string path) {
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


bool FS_check_fasta(std::string& path) {
    std::string line;
    bool valid = false;
    try {
        std::ifstream file(path);
        while (getline(file, line)) {
            if (line.at(0) == '>') valid = true;
        }
        file.close();
        return valid;
    } catch (...) {
        return false;
    }
}


bool FS_create_dir(std::string& path) {
#ifdef USE_BOOST
    return boostFS::create_directories(path);
#endif
}