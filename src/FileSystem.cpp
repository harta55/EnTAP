/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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
#include <sys/stat.h>
#include "config.h"

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
const std::string FileSystem::EXT_TSV  = ".tsv";
const std::string FileSystem::EXT_CSV  = ".csv";
const std::string FileSystem::EXT_LST  = ".lst";
const std::string FileSystem::EXT_CDS  = ".cds";
const std::string FileSystem::EXT_PEP  = ".pep";

const char FileSystem::DELIM_TSV = '\t';
const char FileSystem::DELIM_CSV = ',';
const char FileSystem::FASTA_FLAG = '>';

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


FileSystem::~FileSystem() {
    FS_dprint("Killing object - FileSystem");
    delete_dir(mTempOutpath);
}

FileSystem::FileSystem() {
    // This routine will eventually process entire root directory here and generate
    // hierarchy
    mExeDirectory = "";
    mFinalOutpath = "";
    mTrancriptomeDir = "";
    mTempOutpath = "";

    FS_dprint("Spawn Object - FileSystem");
    set_executable_dir();
    mOriginalWorkingDir = get_cur_dir();
    mRootPath = mOriginalWorkingDir;        // Set root to CWD by default, then change
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
 * Function void print_debug(const std::string &msg)
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
    std::ofstream log_file(mLogFilePath, std::ios::out | std::ios::app);
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
 * @param path          - Path to file to delete
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
    FS_dprint("Deleting file at: " + path);
    return remove(path.c_str()) == 0;
#endif
}


/**
 * ======================================================================
 * Function bool FS_directory_iterate(bool del, std::string &path)
 *
 * Description          - Iterates through a FileSystem directory hierarchy
 *                        performing specified action on files
 *
 * Notes                - WARNING recursive iteration
 *
 * @param iter          - Specify what we would like to do to iterated files
 * @param path          - Absolute path to filesystem directory to iterate through
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
                // If this is a directory, WARNING recursive
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

/**
 * ======================================================================
 * Function std::string set_working_dir(std::string &working_directory)
 *
 * Description          - Sets the current working directory
 *
 * Notes                - POSIX
 *
 * @param working_directory - Working directory we want to set
 *
 * @return              - True (successful), false (unsuccessful)
 *
 * =====================================================================
 */
bool FileSystem::set_working_dir(std::string &working_directory) {
    int ret_code;

    if (working_directory == mCurrentWorkingDir) {
        FS_dprint("Directory is already working: " + working_directory);
        return true;
    }
    if (!file_exists(working_directory)) {
        FS_dprint("Working directory does not exist, creating...");
        if (!create_dir(working_directory)) {
            FS_dprint("ERROR unable to create working directory at: " + working_directory);
            return false;
        }
    }
    FS_dprint("Setting working directory to: " + working_directory);
    ret_code = chdir(working_directory.c_str());

    if (ret_code == 0) {
        mCurrentWorkingDir = working_directory;
        FS_dprint("Success! Working directory set");
        return true;
    } else {
        FS_dprint("ERROR unable to set working directory");
        return false;
    }
}

void FileSystem::set_root_dir(std::string &root) {

    mRootPath     = root;
    mFinalOutpath = PATHS(root, ENTAP_FINAL_OUTPUT);
    mTempOutpath  = PATHS(root, TEMP_DIRECTORY);
    mTrancriptomeDir = PATHS(root, ENTAP_TRANSCRIPTOME_DIR);

    // Make sure directories are created (or already created)
    create_dir(root);
    create_dir(mFinalOutpath);
    delete_dir(mTempOutpath);
    create_dir(mTempOutpath);

    // generate log file
    init_log();
}


/**
 * ======================================================================
 * Function init_log()
 *
 * Description          - Initializes log and debug files
 *
 * Notes                - Sets globals DEBUG_FILE_PATH and LOG_FILE_PATH
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
    debug_file_name = DEBUG_FILENAME + time_date + DEBUG_EXTENSION;

    // Initialize globals
    DEBUG_FILE_PATH = PATHS(mRootPath, debug_file_name);
    mLogFilePath   = PATHS(mRootPath, log_file_name);
    delete_file(DEBUG_FILE_PATH);
    delete_file(mLogFilePath);
    FS_dprint("Start - EnTAP");
}

const std::string &FileSystem::get_root_path() const {
    return mRootPath;
}

/**
 * ======================================================================
 * Function std::string FileSystem::get_file_extension(const std::string &path, bool stripped)
 *
 * Description          - Returns extension of file (path) with or without
 *                        (stripped=TRUE) '.' char
 *                      - Uses Boost FileSystem libraries or Posix
 *
 * Notes                - None
 *
 * @param path          - Absolute path to file
 * @param stripped      - TRUE to remove '.' char from extension, FALSE
 *                        otherwise
 *
 * @return              - Extension of file path
 *
 * =====================================================================
 */
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

/**
 * ======================================================================
 * Function bool FileSystem::copy_file(std::string inpath, std::string outpath, bool overwrite)
 *
 * Description          - Copies file from input path to output path
 *                      - Utilizes Boost FileSystem library or Posix
 *
 * Notes                - None
 *
 * @param inpath        - Absolute path to file we are copying
 * @param outpath       - Absolute path to location where file should be copied
 * @param overwrite     - TRUE to overwrite if outpath file exists already
 *
 * @return              - TRUE if copy was successful
 *
 * =====================================================================
 */
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

    if (!file_exists(inpath)) {
        return false;
    } else if (file_exists(outpath)) {

        if (overwrite) {
            delete_file(outpath);
        } else {
            return false;
        }
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
    STR_REPLACE(path, '\\', '/');
    uint64 dir_pos = path.find_last_of('/');
    uint64 ext_pos = path.find_last_of('.');

    if (dir_pos == std::string::npos && ext_pos == std::string::npos) return path;
    if (ext_pos == std::string::npos) return path.substr(dir_pos + 1);
    if (dir_pos == std::string::npos && with_extension) return path;

    // Filename with first extension (removes trailing extensions)
    if (with_extension) {
        return path.substr(dir_pos + 1);
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
    return this->mFinalOutpath;
}

std::string FileSystem::get_temp_outdir() {
    return this->mTempOutpath;
}

bool FileSystem::download_ftp_file(std::string ftp_path, std::string& out_path) {
    TerminalData terminalData;
    clear_error();

    terminalData.print_files = false;
    terminalData.suppress_std_err = false;
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
    clear_error();

    terminalData.print_files = false;
    terminalData.suppress_std_err = false;
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
        case ENT_FILE_TAR_GZ:
            terminal_cmd =
                "tar -xzf " + in_path + " -C " + out_dir;
            break;

        case ENT_FILE_GZ:
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

/**
 * ======================================================================
 * Function bool FileSystem::rename_file(std::string &in, std::string &out)
 *
 * Description          - Renames input file using Boost FileSystem library
 *                        or POSIX
 *                      - Will move files to out path as well
 *
 * Notes                - None
 *
 * @param in           - Absolute path to file we are renaming
 * @param out          - Absolute path to renamed file (will be moved if in different
 *                       directory)
 *
 * @return              - TRUE if copy was successful
 *
 * =====================================================================
 */
bool FileSystem::rename_file(std::string &in, std::string &out) {
    FS_dprint("Moving/renaming file: " + in );
    if (!file_exists(in)) {
        FS_dprint("ERROR: rename failed, file does not exist");
        return false;
    }
    try {
#ifdef USE_BOOST
        boostFS::rename(in, out);
        FS_dprint("Success!");
        return true;
#else
        rename(in.c_str(), out.c_str());
#endif
    } catch (...) {
        FS_dprint("ERROR: rename failed!");
        return false;
    }
    return true;
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
    mErrorMsg = err_msg;
}

std::string FileSystem::get_error() {
    return "\n" + mErrorMsg;
}

void FileSystem::format_stat_stream(std::stringstream &stream, std::string title) {
    stream<<std::fixed<<std::setprecision(2);
    stream << SOFTWARE_BREAK << title << '\n' << SOFTWARE_BREAK;

}

std::string FileSystem::get_extension(FileSystem::ENT_FILE_TYPES type) {
    switch (type) {

        case ENT_FILE_DELIM_TSV:
        case ENT_FILE_GENE_ENRICH_GO_TERM:
        case ENT_FILE_GENE_ENRICH_EFF_LEN:
            return EXT_TSV;

        case ENT_FILE_DELIM_CSV:
            return EXT_CSV;

        case ENT_FILE_XML:
            return EXT_XML;

        case ENT_FILE_FASTA_FNN:
            return EXT_FNN;

        case ENT_FILE_FASTA_FAA:
            return EXT_FAA;

        default:
            FS_dprint("ERROR unhandled extension type: " + std::to_string(type));
            return "";
    }
}

bool FileSystem::create_transcriptome_dir() {
    return create_dir(mTrancriptomeDir);
}

std::string FileSystem::get_trancriptome_dir() {
    return mTrancriptomeDir;
}

/**
 * ======================================================================
 * Function std::string FileSystem::set_executable_dir()
 *
 * Description          - Gets the directory that the EnTAP executable was detected
 *                        in
 *
 * Notes                - Only implemented for UNIX systems now!
 *
 *
 * @return              - Path to directory that EnTAP exe was found
 * ======================================================================
 */
void FileSystem::set_executable_dir() {
    mExeDirectory = get_exe_dir();
}

std::string FileSystem::get_exe_dir() {
    try {
        char buff[1024];
        ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
        if (len != -1) {
            buff[len] = '\0';
            std::string path = std::string(buff);
            return path.substr(0,path.find_last_of('/'));
        }
    } catch (...) {
        return "";
    }
    return "";
}

bool FileSystem::is_url(std::string &url) {
    return ((url.find("http") != std::string::npos && (url.find("www") != std::string::npos)));
}

void FileSystem::clear_error() {
    mErrorMsg = "";
}

/**
 * ======================================================================
 * Function bool FileSystem::move_dir(std::string &original, std::string &final)
 *
 * Description          - Moves input file using Boost FileSystem library
 *                        or POSIX (calls rename_file())
 *
 * Notes                - None
 *
 * @param original     - Absolute path to dir we want to move
 * @param final        - Absolute path to where we are moving the directory
 *
 * @return              - TRUE if copy was successful
 *
 * =====================================================================
 */
bool FileSystem::move_dir(std::string &original, std::string &final) {
    return rename_file(original, final);    // rename_file checks for directory existence
}

/**
* ======================================================================
* Function std::string ModInterpro::format_for_csv_parser(const std::string &input_path, std::string &output_path,
*                                                      uint16 col_num)
*
* Description          - Formats tsv file to be read easier with current
*                        lib (complains if tabs are not perfect)
*                      - This will add in tabs to keep everything consistent
*                     - Temporary...
*
*
* Notes                - None
*
*
* @param input_path    - Absolute path to file we want to format
* @param output_path   - Path to temporary formatted file
* @param col_num       - Number of columns in input/output
*
* @return              - FALSE if error occurred, TRUE if sucesful
*
* =====================================================================
*/
bool FileSystem::format_for_csv_parser(const std::string &input_path, std::string &output_path, uint16 col_num) {
    std::string path_temp;
    std::string line;
    std::string input_filename;
    uint16      tab_ct;
    bool        ret = true;

    path_temp = input_path;

    input_filename = get_filename(path_temp, false) + "_temp";
    path_temp      = PATHS(mTempOutpath, input_filename);
    delete_file(path_temp);

    try {
        std::ifstream file_in(input_path);
        std::ofstream file_temp(path_temp, std::ios::out | std::ios::app);
        while(std::getline(file_in, line)) {
            if (line.empty()) continue;
            file_temp << line;
            tab_ct = (uint16)std::count(line.begin(), line.end(), '\t');
            while (tab_ct < col_num - 1) {
                file_temp<<'\t';
                tab_ct++;
            }
            file_temp << std::endl;
        }
        file_in.close();
        file_temp.close();
    } catch (std::exception &e) {
        FS_dprint("ERROR format_for_csv_parser: " + std::string(e.what()));
        ret = false;
    }


    output_path = path_temp;        // Set output
    return ret;
}

const std::string &FileSystem::getMOriginalWorkingDir() const {
    return mOriginalWorkingDir;
}
