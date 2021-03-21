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

#ifndef ENTAP_FILESYSTEM_H
#define ENTAP_FILESYSTEM_H


//*********************** Includes *****************************
#include "common.h"
#include "config.h"
#include "TerminalCommands.h"
#include "EntapGlobals.h"
#ifdef USE_BOOST
#include <boost/filesystem.hpp>
#endif
//**************************************************************

//******************** Defines/Macros **************************
#ifdef USE_BOOST
#define PATHS(x,y)      (boostFS::path(x) / boostFS::path(y)).string()
#else
#define PATHS(x,y)      ((x) + "/" + (y))
#endif
//**************************************************************


//***************** Global Prototype Functions *****************
void FS_dprint(const std::string&);
//**************************************************************

class FileSystem {

public:

    typedef enum {
        ENT_FILE_UNUSED=0,
        ENT_FILE_DELIM_TSV,
        ENT_FILE_DELIM_CSV,
        ENT_FILE_FASTA_FAA,
        ENT_FILE_FASTA_FNN,
        ENT_FILE_GENE_ENRICH_EFF_LEN,           // Gene enrichment format for gene ontology
                                                // gene ID and effective length, tab delim
        ENT_FILE_GENE_ENRICH_GO_TERM,           // Gene ID and Go terms, tab delim
                                                // new row for every go term for each gene
        ENT_FILE_OUTPUT_FORMAT_MAX,     // File types above this are supported for data output

        ENT_FILE_XML,                   // Not yet supported for output format
        ENT_FILE_TAR_GZ,
        ENT_FILE_GZ,
        ENT_FILE_ZIP,

        ENT_FILE_MAX

    } ENT_FILE_TYPES;

    typedef enum {
        FILE_STATUS_UNUSED     = (1 << 0),
        FILE_STATUS_EMPTY      = (1 << 1),      // File is empty
        FILE_STATUS_READ_ERR   = (1 << 2),      // Error reading file
        FILE_STATUS_PATH_ERR   = (1 << 3),      // FIle does not exist

        FILE_STATUS_MAX        = (1 << 15)

    } ENT_FILE_STATUS;

    const std::string &getMOriginalWorkingDir() const;

    typedef enum {

        FILE_ITER_DELETE=0,
        FILE_ITER_DELETE_EMPTY,
        FILE_ITER_PRINT

    } ENT_FILE_ITER;


    FileSystem();
    ~FileSystem();
    void set_root_dir(std::string &root_dir);
    bool set_working_dir(std::string &working_dir);
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
    bool move_dir(std::string& original, std::string &final);
    void delete_dir(std::string&);

    bool create_transcriptome_dir();

    const std::string &get_root_path() const;
    std::string get_file_extension(const std::string&, bool);
    std::string get_filename(std::string&, bool);
    static std::string get_cur_dir();
    static std::string get_exe_dir();
    std::vector<std::string> list_to_vect(char, std::string&);
    std::string get_final_outdir();
    std::string get_temp_outdir();
    bool rename_file(std::string& in, std::string& out);
    uint16 get_file_status(std::string &path);
    std::string print_file_status(uint16 status,std::string& path);
    std::string get_error();
    std::string get_extension(ENT_FILE_TYPES type);
    std::string get_trancriptome_dir();
    bool format_for_csv_parser(const std::string &input_path, std::string &output_path, uint16 col_num);

    bool download_ftp_file(std::string,std::string&);
    bool decompress_file(std::string &in_path, std::string &out_dir, ENT_FILE_TYPES);
    bool is_url(std::string &url);

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
    static const std::string EXT_TSV;
    static const std::string EXT_CSV;
    static const std::string EXT_LST;
    static const std::string EXT_PEP;
    static const std::string EXT_CDS;

    static const char        DELIM_TSV;
    static const char        DELIM_CSV;
    static const char        FASTA_FLAG;
#ifndef UNIT_TESTS
private:
#endif
    //****************** Private Functions *********************
    void init_log();
    void set_error(std::string err_msg);
    void set_executable_dir();
    void clear_error();
    //**********************************************************

    //**************** Private Const Variables *****************
    const std::string LOG_FILENAME              = "log_file"; // Filename for EnTAP log file (statistics files)
    const std::string LOG_EXTENSION             = EXT_TXT; // Extension for EnTAP statistics file
    const std::string DEBUG_EXTENSION           = EXT_TXT; // Extension for EnTAP debug file
#ifdef UNIT_TESTS
    const std::string DEBUG_FILENAME            = "unit_testing";
#else
    const std::string DEBUG_FILENAME            = "debug"; // Filename for EnTAP debug file
#endif
    const std::string ENTAP_FINAL_OUTPUT        = "final_results/"; // Directory name for final output annotations directory
    const std::string ENTAP_TRANSCRIPTOME_DIR   = "transcriptomes/"; // Directory name for transcriptome directory (frame selected, expression analysis)
    const std::string TEMP_DIRECTORY            = "temp/"; // Directory name for 'temp' directory  (deleted once EnTAP exits)

    // Log file specifics
    const std::string SOFTWARE_BREAK            = "------------------------------------------------------\n";
    //**********************************************************

    //****************** Private Variables *********************
    std::string mLogFilePath;
    std::string mRootPath;        // Root EnTAP output directory
    std::string mFinalOutpath;    // Path to final files after entap has finished
    std::string mTrancriptomeDir; // Absolute path to EnTAP transcriptome directory
    std::string mTempOutpath;     // Temp directory for EnTAP usage
    std::string mErrorMsg;        // String containing error mMessage from execution
    std::string mExeDirectory;    // Directory of the EnTAP executable
    std::string mOriginalWorkingDir;     // Original working directory;
    std::string mCurrentWorkingDir;      // Current working directory
    //**********************************************************
};



#endif //ENTAP_FILESYSTEM_H
