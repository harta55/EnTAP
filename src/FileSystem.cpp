/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
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
#include <sys/stat.h>
#include "config.h"
#include "TerminalCommands.h"

#ifdef USE_BOOST
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/time_clock.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#else
#include <dirent.h> // POSIX directory iterator
#include <zconf.h>
#endif

const std::string FileSystem::EXT_TXT  = ".txt";
const std::string FileSystem::EXT_ERR  = ".err";
const std::string FileSystem::EXT_OUT  = ".out";
const std::string FileSystem::EXT_BAM  = ".bam";
const std::string FileSystem::EXT_SAM  = ".sam";
const std::string FileSystem::EXT_FAA  = ".faa";
const std::string FileSystem::EXT_FNN  = ".fnn";
const std::string FileSystem::EXT_XML  = ".xml";
const std::string FileSystem::EXT_DMND = ".dmnd";
const std::string FileSystem::EXT_STD  = "_std";

const char FileSystem::DELIM_TSV= '\t';

// Removed for older compilers, may bring back
#if 0
void FileSystem::open_out(std::string &path, std::ofstream &ofstream) {
    ofstream = std::ofstream(path,std::ios::out | std::ios::app);
    if (!ofstream.is_open()) {
        throw ExceptionHandler("Error opening file: " + path,
             ERR_ENTAP_FILE_IO);
    }
}
#endif


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
 *                      - Adds bo to each entry
 *
 * Notes                - None
 *
 * @param msg           - Message to be sent to debug file
 * @return              - None
 *
 * =====================================================================
 */
void FS_dprint(const std::string &msg) {
    std::ofstream debug_file(DEBUG_FILE_PATH, std::ios::out | std::ios::app);

    debug_file << get_cur_time() << ": " + msg << std::endl;
    debug_file.close();
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
//    return ( access( path.c_str(), F_OK ) != -1 );
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
    return remove(path.c_str()) == 0;
#endif
}


/**
 * ======================================================================
 * Function bool FS_directory_iterate(bool del, std::string &path)
 *
 * Description          - WARNING recursive iteration
 *
 * Notes                - None
 *
 * @param path          - Path to directory
 * @param del           - True/false to delete files that are empty (unused)
 *
 * @return              - True/false if successful
 * ======================================================================
 */
bool FileSystem::directory_iterate(ENT_FILE_ITER iter, std::string &path) {
//    FS_dprint("Iterating through directory: " + path);
#ifdef USE_BOOST
    if (!file_exists(path)) return false;
    try {
        for (boostFS::recursive_directory_iterator it(path), end; it != end; ++it) {
            if (!boostFS::is_directory(it->path())) {
                // Is file
                if (file_empty(it->path().string()) && iter == FILE_ITER_DELETE_EMPTY) {
                    delete_file(it->path().string());
                    FS_dprint("Deleted: " + it->path().string());
                }
            }
        }
    } catch (...) {
        return false;
    }
    return true;
#else   // POSIX
    struct dirent *entry;
    DIR *dp;
    std::string file_path;
    try {

        dp = opendir(path.c_str());
        if (dp == NULL) {
            FS_dprint("opendir: Path could not be read: " + path);
            return false;
        }

        while ((entry = readdir(dp)) != NULL) {
            if (entry->d_name[0] == '.') continue;
            file_path = PATHS(path, std::string(entry->d_name));
            if (entry->d_type == DT_DIR) {
                // Directory, WARNING recursive
                directory_iterate(FILE_ITER_DELETE, file_path);     // Recursive call
            } else {
                // This is a file decide what we want to do
                switch (iter) {
                    case FILE_ITER_DELETE:
                        delete_file(file_path);
                        break;
                    case FILE_ITER_DELETE_EMPTY:
                        if (file_empty(file_path)) delete_file(file_path);
                        break;
                    default:
                        break;
                }
            }
        }
        closedir(dp);
        remove(path.c_str());
        return true;
    } catch (...) {return false;};
#endif
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
#else
    return mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != -1;
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
#else   // POSIX
    directory_iterate(FILE_ITER_DELETE, path);
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
#else // POSIX
    char cwd[PATH_MAX];
    getcwd(cwd, sizeof(cwd));
    return std::string(cwd);
#endif

}

FileSystem::~FileSystem() {
    FS_dprint("Killing object - FileSystem");
    delete_dir(_temp_outpath);
}

FileSystem::FileSystem(std::string &root) {
    // This routine will process entire root directory here and generate
    // hierarchy
    FS_dprint("Spawn Object - FileSystem");

    _root_path     = root;
    _final_outpath = PATHS(root, ENTAP_FINAL_OUTPUT);
    _temp_outpath  = PATHS(root, TEMP_DIRECTORY);

    // Make sure directories are created (or already created)
    create_dir(root);
    create_dir(_final_outpath);
    create_dir(_temp_outpath);

    // generate log file
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

    // Will not work on older GCC compilers with a known bug in put_time
    /*
    std::chrono::time_point<std::chrono::system_clock> now;
    std::time_t                                        now_time;
    std::tm                                            now_tm;
    now      = std::chrono::system_clock::now();
    now_time = std::chrono::system_clock::to_time_t(now);
    now_tm   = *std::localtime(&now_time);
    ss << std::put_time(&now_tm, "_%YY%mM%dD-%Hh%Mm%Ss");
     */

    time_t theTime = time(nullptr);
    struct tm *aTime = localtime(&theTime);

    ss << "_" << aTime->tm_year+1900 << "Y"
              << aTime->tm_mon+1  << "M"
              << aTime->tm_mday << "D-"
              << aTime->tm_hour << "h"
              << aTime->tm_min  << "m"
              << aTime->tm_sec  << "s";
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

// Stipped without '.'
std::string FileSystem::get_file_extension(const std::string &path, bool stripped) {
#ifdef USE_BOOST
    boostFS::path bpath(path);
    if (stripped) { // remove period
        return bpath.extension().string().substr(1);
    } else {
        return bpath.extension().string();
    }
#else
    if (path.find_last_of('.') != std::string::npos) {
        if (stripped) {
            return path.substr(path.find_last_of('.') + 1);
        } else {
            return path.substr(path.find_last_of('.'));
        }
    } else {
        return "";
    }
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
#else
    if (file_exists(outpath) && !overwrite) {
        return false;
    } else if (file_exists(outpath)) {
        delete_file(outpath);
    }
    std::ifstream  in(inpath, std::ios::binary);
    std::ofstream  out(outpath,   std::ios::binary);
    out << in.rdbuf();
    in.close();
    out.close();
    return true;
#endif
}

std::string FileSystem::get_filename(std::string &path, bool with_extension) {
    if (path.empty()) return "";
#ifdef USE_BOOST
    if(with_extension) {
        boostFS::path boost_path(path);
        return boost_path.filename().string();
    } else {
        std::string file_name;
        boostFS::path path(file_name);
        while (path.has_extension()) path = path.stem();
        file_name = path.string();
        return file_name;
    }
#else
    unsigned long dir_pos = path.find_last_of('/');
    unsigned long ext_pos = path.find_last_of('.');

    if (dir_pos == std::string::npos && ext_pos == std::string::npos) return path;
    if (ext_pos == std::string::npos) return path.substr(dir_pos + 1);
    if (dir_pos == std::string::npos && with_extension) return path;

    // Filename with first extension (removes trailing extensions)
    if (with_extension) {
        return path.substr(dir_pos + 1, ext_pos - dir_pos - 1);
    } else {
        // Just get stem
        if (dir_pos == std::string::npos) {
            return path.substr(0, path.find('.'));
        } else {
            return path.substr(dir_pos+1, path.find('.', dir_pos) - dir_pos - 1);
        }
    }
#endif
}


std::string FileSystem::get_final_outdir() {
    return this->_final_outpath;
}

std::string FileSystem::get_temp_outdir() {
    return this->_temp_outpath;
}

bool FileSystem::download_ftp_file(std::string ftp_path, std::string& out_path) {
    TerminalData terminalData;

    terminalData.print_files = false;
    terminalData.base_std_path = "";


    FS_dprint("Downloading FTP file at: " + ftp_path);

#ifdef USE_CURL
    FS_dprint("Using CURL...");
    CURL *curl;
    CURLcode res;
    FILE *out_file;

    curl = curl_easy_init();
    if (curl) {
        out_file = fopen(out_path.c_str(),"wb");
        curl_easy_setopt(curl, CURLOPT_URL, ftp_path.c_str());
        // We dont need a callback for now
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, NULL);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, out_file);
        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
        fclose(out_file);

        if(CURLE_OK != res) {
            // Failed
            FS_dprint("CURL download has failed");
            return false;
        }
    } else {
        FS_dprint("CURL has failed!");
        return false;
    }
    return true;
#else
    // Not compiled with CURL usage, use terminal command
    FS_dprint("Using wget terminal command...");

    std::string terminal_cmd =
        "wget -O " + out_path + " " + ftp_path;

    terminalData.command = terminal_cmd;

    if (TC_execute_cmd(terminalData) == 0) {
        FS_dprint("Success, file saved to: " + out_path);
        return true;
    } else {
        set_error("Unable to get file\n" + terminalData.err_stream);
        return false;
    }
#endif
}

bool FileSystem::decompress_file(std::string &in_path, std::string &out_dir, ENT_FILE_TYPES type) {
    TerminalData terminalData;

    terminalData.print_files = false;
    terminalData.base_std_path = "";

    FS_dprint("Decompressing file at: " + in_path);
    if (!file_exists(in_path)) {
        FS_dprint("File does not exist!");
        return false;
    }
#ifdef USE_ZLIB
    FS_dprint("Using ZLIB...");
    FS_dprint("ZLIP NOT SUPPORTED");
    return false;
#else
    // Not compiled with ZLIB usage, use terminal command
    FS_dprint("Using terminal command...");
    std::string terminal_cmd;

    switch (type) {
        case FILE_TAR_GZ:
            terminal_cmd =
                "tar -xzf " + in_path + " -C " + out_dir;
            break;

        case FILE_GZ:
            terminal_cmd =
                "gunzip -c " + in_path + " > " + out_dir; //outdir will be outpath in this case
            break;
        default:
            return false;
    }

    terminalData.command = terminal_cmd;

    if (TC_execute_cmd(terminalData) == 0) {
        FS_dprint("Success! Exported to: " + out_dir);
        return true;
    } else {
        set_error("Unable to decompress file\n" + terminalData.err_stream);
        return false;
    }
#endif
}

bool FileSystem::rename_file(std::string &in, std::string &out) {
    FS_dprint("Moving/renaming file: " + in );
    if (!file_exists(in)) return false;
    try {
#ifdef USE_BOOST
        boostFS::rename(in, out);
        FS_dprint("Success!");
        return true;
#else
        rename(in.c_str(), out.c_str());
#endif
    } catch (...) {
        FS_dprint("rename failed!");
        return false;
    }
    return false;
}

uint16 FileSystem::get_file_status(std::string &path) {
    uint16 file_status = 0;

    if (!file_exists(path)) {
        file_status |= FILE_STATUS_PATH_ERR;
    }

    if (file_empty(path)) {
        file_status |= FILE_STATUS_EMPTY;
    }

    if (!file_test_open(path)) {
        file_status |= FILE_STATUS_READ_ERR;
    }

    return file_status;
}

std::string FileSystem::print_file_status(uint16 status, std::string& path) {
    std::string err_msg;

    err_msg = "\nFile Status at: " + path;
    if (status == 0) return err_msg += "\n    No Error";

    // Check if file empty
    if ((status & FILE_STATUS_EMPTY) != 0) {
        err_msg += "\n    Error: File empty";
    }

    // Check if file path is invalid
    if ((status & FILE_STATUS_PATH_ERR) != 0) {
        err_msg += "\n    Error: File could not be found";
    }

    // Check if file could not be opened
    if ((status & FILE_STATUS_READ_ERR) != 0) {
        err_msg += "\n    Error: File could not be opened";
    }

    return err_msg;
}

void FileSystem::set_error(std::string err_msg) {
//    FS_dprint(err_msg);
    _err_msg = err_msg;
}

std::string FileSystem::get_error() {
    return "\n" + _err_msg;
}

bool FileSystem::print_headers(std::ofstream& file_stream, std::vector<ENTAP_HEADERS> &headers, char delim) {
    for (ENTAP_HEADERS &header: headers) {
        if (ENTAP_HEADER_INFO[header].print_header) {
            file_stream << ENTAP_HEADER_INFO[header].title << delim;
        }
    }
    file_stream << std::endl;
    return true;
}

void FileSystem::format_stat_stream(std::stringstream &stream, std::string title) {
    stream<<std::fixed<<std::setprecision(2);
    stream << SOFTWARE_BREAK << title << '\n' << SOFTWARE_BREAK;

}
