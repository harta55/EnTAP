/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
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
//**************************************************************


// Keeping global for now
void FS_dprint(const std::string&);

typedef enum {
    FILE_TAR_GZ,
    FILE_ZIP

} ENT_FILE_TYPES;

//***************** Global Prototype Functions *****************
class FileSystem {

public:

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
    bool directory_iterate(bool, std::string&);
    bool check_fasta(std::string&);
    bool create_dir(std::string&);
    void delete_dir(std::string&);
    const std::string &get_root_path() const;
    std::string get_file_extension(const std::string&, bool);
    void get_filename_no_extensions(std::string &);
    std::string get_filename(std::string);
    static std::string get_cur_dir();
    std::vector<std::string> list_to_vect(char, std::string&);
    std::string get_final_outdir();
    std::string get_temp_outdir();

    bool download_ftp_file(std::string&,std::string&);
    bool extract_file(std::string &in_path, std::string &out_path);

//**************************************************************
    static const std::string EXT_TXT ;
    static const std::string EXT_ERR ;
    static const std::string EXT_OUT ;
    static const std::string EXT_BAM ;
    static const std::string EXT_FAA ;
    static const std::string EXT_FNN ;
    static const std::string EXT_DMND;
    static const std::string EXT_XML;

private:
    void init_log();

    const std::string LOG_FILENAME   = "log_file";
    const std::string LOG_EXTENSION  = ".txt";
    const std::string DEBUG_FILENAME = "debug";
    const std::string ENTAP_FINAL_OUTPUT    = "final_results/";
    const std::string TEMP_DIRECTORY        = "temp/";
    std::string _root_path;     // Root EnTAP output directory
    std::string _final_outpath; // Path to final files after entap has finished
    std::string _temp_outpath;  // Temp directory for EnTAP usage
};



#endif //ENTAP_FILESYSTEM_H
