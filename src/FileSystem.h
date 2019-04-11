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

#ifndef ENTAP_FILESYSTEM_H
#define ENTAP_FILESYSTEM_H


//*********************** Includes *****************************
#include "common.h"
#include "TerminalCommands.h"
#include "EntapGlobals.h"
//**************************************************************


// Keeping global for now
void FS_dprint(const std::string&);

//***************** Global Prototype Functions *****************
class FileSystem {

public:

    typedef enum {
        FILE_TAR_GZ,
        FILE_GZ,
        FILE_ZIP,
        FILE_DELIMINATED,
        FILE_FASTA,

        FILE_MAX

    } ENT_FILE_TYPES;

    typedef enum {
        FILE_STATUS_UNUSED     = (1 << 0),
        FILE_STATUS_EMPTY      = (1 << 1),      // File is empty
        FILE_STATUS_READ_ERR   = (1 << 2),      // Error reading file
        FILE_STATUS_PATH_ERR   = (1 << 3),      // FIle does not exist

        FILE_STATUS_MAX        = (1 << 15)

    } ENT_FILE_STATUS;

    typedef enum {

        FILE_ITER_DELETE=0,
        FILE_ITER_DELETE_EMPTY,
        FILE_ITER_PRINT

    } ENT_FILE_ITER;

    struct CreateFileCmd {

        std::stringstream *stringstream;
        std::string        file_path;
        std::vector<ENTAP_HEADERS> headers;
        ENT_FILE_TYPES     type;
        char               delim;
    };

    FileSystem(std::string&);
    ~FileSystem();
    //void open_out(std::string &, std::ofstream &);
    bool file_is_open(std::ofstream&);
    void close_file(std::ofstream&);
    void print_stats(std::string &msg);
    bool file_test_open(std::string&);
    bool file_exists(std::string);
    bool file_empty(std::string);
    bool file_no_lines(std::string);
    bool delete_file(std::string);
    bool copy_file(std::string, std::string, bool);
    bool directory_iterate(ENT_FILE_ITER, std::string&);
    bool check_fasta(std::string&);
    bool create_dir(std::string&);
    void delete_dir(std::string&);
    const std::string &get_root_path() const;
    std::string get_file_extension(const std::string&, bool);
    std::string get_filename(std::string&, bool);
    static std::string get_cur_dir();
    std::vector<std::string> list_to_vect(char, std::string&);
    std::string get_final_outdir();
    std::string get_temp_outdir();
    bool rename_file(std::string& in, std::string& out);
    uint16 get_file_status(std::string &path);
    std::string print_file_status(uint16 status,std::string& path);
    std::string get_error(void);

    bool download_ftp_file(std::string,std::string&);
    bool decompress_file(std::string &in_path, std::string &out_dir, ENT_FILE_TYPES);

    bool print_headers(std::ofstream &file_stream, std::vector<ENTAP_HEADERS> &headers, char delim);
    void format_stat_stream(std::stringstream &stream, std::string title);

//**************************************************************
    static const std::string EXT_TXT ;
    static const std::string EXT_ERR ;
    static const std::string EXT_OUT ;
    static const std::string EXT_BAM ;
    static const std::string EXT_SAM ;
    static const std::string EXT_FAA ;
    static const std::string EXT_FNN ;
    static const std::string EXT_DMND;
    static const std::string EXT_XML;
    static const std::string EXT_STD;

    static const char        DELIM_TSV;

private:
    void init_log();
    void set_error(std::string err_msg);

    const std::string LOG_FILENAME   = "log_file";
    const std::string LOG_EXTENSION  = ".txt";
    const std::string DEBUG_FILENAME = "debug";
    const std::string ENTAP_FINAL_OUTPUT    = "final_results/";
    const std::string TEMP_DIRECTORY        = "temp/";
    const std::string SOFTWARE_BREAK = "------------------------------------------------------\n";

    std::string _root_path;     // Root EnTAP output directory
    std::string _final_outpath; // Path to final files after entap has finished
    std::string _temp_outpath;  // Temp directory for EnTAP usage
    std::string _err_msg;
};



#endif //ENTAP_FILESYSTEM_H
