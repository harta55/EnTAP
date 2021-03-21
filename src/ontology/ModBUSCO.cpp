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

#include <csv.h>
#include "ModBUSCO.h"
#include "../QuerySequence.h"
#include "../QueryData.h"
#include "../ExceptionHandler.h"

std::vector<ENTAP_HEADERS> ModBUSCO::DEFAULT_HEADERS = {
    ENTAP_HEADER_ONT_BUSCO_ID,
    ENTAP_HEADER_ONT_BUSCO_STATUS,
    ENTAP_HEADER_ONT_BUSCO_LENGTH,
    ENTAP_HEADER_ONT_BUSCO_SCORE
};

ModBUSCO::ModBUSCO(std::string &in_hits, std::string &ont_out, EntapDataPtrs &entap_data) :
        AbstractOntology(in_hits, ont_out, entap_data, "BUSCO", DEFAULT_HEADERS) {
    FS_dprint("Spawn Object - ModBUSCO");

    // Pull inputs from user
    mExePath       = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_BUSCO_EXE);
    mBuscoDatabase = mpUserInput->get_user_input<ent_input_str_t>(INPUT_FLAG_BUSCO_DATABASE);
    mEval          = mpUserInput->get_user_input<ent_input_fp_t >(INPUT_FLAG_BUSCO_EVAL);
    mBlastp ? mRunType = BUSCO_RUN_TYPE_PROT : mRunType= BUSCO_RUN_TYPE_TRAN;
    mBuscoVersion  = BUSCO_VERSION_UNKNOWN;

    mOutputDirectoryTag = "busco_" + mpFileSystem->get_filename(mBuscoDatabase, false);
    mOutputRunDir    = PATHS(FileSystem::get_cur_dir(), mOutputDirectoryTag);

    set_version();  // Set BUSCO software version
    mFinalTablePath = get_final_table_path(mOutputRunDir, mBuscoDatabase, mBuscoVersion);
}

ModBUSCO::~ModBUSCO() {
    FS_dprint("Killing object - ModBUSCO");
}

EntapModule::ModVerifyData ModBUSCO::verify_files() {
    ModVerifyData ret  = ModVerifyData();

    FS_dprint("Looking for BUSCO file at: " + mFinalTablePath);

    // Is our final table present?
    if (!mpFileSystem->file_exists(mFinalTablePath) || mpFileSystem->file_empty(mFinalTablePath)) {
        ret.files_exist = false;
        FS_dprint("BUSCO output file not found, continuing with execution...");
    } else {
        // Our file exists and is not empty
        FS_dprint("BUSCO file found! Skipping execution");
        ret.output_paths.push_back(mFinalTablePath);
        ret.files_exist = true;
    }

    return ret;
}

bool ModBUSCO::is_executable(std::string &exe) {
    std::string test_command;
    TerminalData terminalData;

    test_command = exe + " --help";

    terminalData.command = test_command;
    terminalData.print_files = false;
    terminalData.suppress_std_err = false;

    return TC_execute_cmd(terminalData) == 0;
}

void ModBUSCO::execute() {
    std::string  database_out_path;     // Base path that BUSCO will output files to
                                        // WARNING!!! BUSCO v3 prepends "run_" to this
    std::string  std_out;
    std::string  cmd;
    TerminalData terminalData;

    // first we want to set the version
    set_version();

    // Is version supported?
    if (!is_version_supported(mBuscoVersion)) {
        // NO version not supported by EnTAP
        throw ExceptionHandler("ERROR unsupported BUSCO version: " + print_version(), ERR_ENTAP_VERSION_UNSUPPORTED_BUSCO);
    } else {
        // YES Version is supported
        // Execute BUSCO analysis against database input by user
        // WARNING database will be  downloaded if not already locally installed
        // WARNING BUSCO v3 and v4 only output to the current working directory

        switch (mBuscoVersion) {

            case BUSCO_VERSION_3:
                break;

            case BUSCO_VERSION_4:
                FS_dprint("Executing BUSCO for version 4");
                // WARNING version 4 has separate ini file for E-val
                std_out = mpFileSystem->get_filename(mBuscoDatabase, false) + "_" + FileSystem::EXT_STD;
                database_out_path = mpFileSystem->get_filename(mBuscoDatabase, false);

                cmd =
                        mExePath + " " +
                        BUSCO_INPUT_IN       + " " + mInputTranscriptome + " " + /* Input Transcriptome*/
                        BUSCO_INPUT_OUTPUT   + " " + database_out_path   + " " + /* Output base path */
                        BUSCO_INPUT_RUN_TYPE + " " + mRunType            + " " + /* Run type (tran or prot) */
                        BUSCO_INPUT_DATABASE + " " + mBuscoDatabase      + " " + /* Database */
                        BUSCO_INPUT_CPU      + " " + std::to_string(mThreads) + " " +
                        BUSCO_INPUT_EVAL     + " " + std::to_string(mEval);

                if (TC_execute_cmd(terminalData) != 0) {
                    // ERROR in run, cleanup
                    mpFileSystem->delete_dir(mOutputRunDir);
                    FS_dprint("BUSCO STD OUT:\n" + terminalData.out_stream);
                    throw ExceptionHandler("Error in running BUSCO against database at: " +
                                           mBuscoDatabase + "\nBUSCO Error:\n" + terminalData.err_stream, ERR_ENTAP_RUN_BUSCO);                }

                break;

            default:
                break;
        }
    }
}

void ModBUSCO::parse() {
    std::stringstream stats_stream;
    std::string temp_reformatted;       // Temporary path to reformatted tsv file
    uint16 file_status=0;

    FS_dprint("Beginning to parse BUSCO file at: " + mFinalTablePath);

    // Ensure file is valid
    file_status = mpFileSystem->get_file_status(mFinalTablePath);
    if (file_status != 0) {
        throw ExceptionHandler(mpFileSystem->print_file_status(file_status,mFinalTablePath),
                               ERR_ENTAP_PARSE_BUSCO);
    }

    // Need to reformat our TSV file since our parser does not like 'clean' tsv files...temporary need to change
    if (!mpFileSystem->format_for_csv_parser(mFinalTablePath, temp_reformatted, BUSCO_COLUMN_NUM)) {
        // ERROR unable to reformat file
        throw ExceptionHandler("ERROR unable to reformat BUSCO file", ERR_ENTAP_PARSE_BUSCO);
    }

    // File valid, continue
    mpFileSystem->format_stat_stream(stats_stream, "Annotation Completeness - BUSCO");

    // Parse through library
    std::string busco_id;
    std::string seq_status;
    std::string sequence_id;
    fp64 busco_score;
    uint32 seq_length;
    vect_str_t missing_buscos;
    QuerySequence::BuscoResults buscoResults;
    QuerySequence  *querySequence;

    io::CSVReader<BUSCO_COLUMN_NUM, io::trim_chars<' '>, io::no_quote_escape<'\t'>, io::throw_on_overflow,io::single_and_empty_line_comment<'#'>>
            in(temp_reformatted);
    while (in.read_row(busco_id, seq_status, sequence_id, busco_score, seq_length)) {

        buscoResults = {};

        // Check if this is a missing busco
        if (sequence_id.empty()) {

            // YES, record
            missing_buscos.push_back(busco_id);

        } else {
            // NO, Check if we recognize this sequence ID
            querySequence = mpQueryData->get_sequence(sequence_id);
            if (querySequence == nullptr) {
                throw ExceptionHandler("ERROR unable to find BUSCO sequence: " + sequence_id + " in user transcriptome",
                                       ERR_ENTAP_PARSE_BUSCO);
            }

            // YES, sequence good. record data
            buscoResults.status = seq_status;
            buscoResults.length = seq_length;
            buscoResults.length_str = std::to_string(seq_length);
            buscoResults.score  = busco_score;
            buscoResults.score_str = float_to_string(busco_score);
            buscoResults.busco_id  = busco_id;

            querySequence->add_alignment(GENE_ONTOLOGY, ONT_BUSCO, buscoResults, mBuscoDatabase);

        }
    } // END WHILE

    FS_dprint("Success! Parsing complete");
}

bool ModBUSCO::set_version() {
    std::string test_command;
    TerminalData terminalData;
    int err_code;                  // Error code from terminal command
    bool ret=true;                 // Able to pull version info

    FS_dprint("Setting BUSCO version...");

    try {
        test_command = mExePath + " --version";

        terminalData.command = test_command;
        terminalData.print_files = true;
        terminalData.suppress_std_err = false;

        err_code = TC_execute_cmd(terminalData);

        // Check if we were able to successfully execute command
        if (err_code == 0) {
            // Able to execute command, parse version from output
            FS_dprint("BUSCO raw version:" + terminalData.out_stream + " parsing...");

            trim(terminalData.out_stream);  // trim spaces

            // Should be printed to terminal like 'BUSCO 4.0.2'
            std::string token = "BUSCO ";   // MUST at beginning of string
            if (terminalData.out_stream.find(token) == 0) {
                auto pos1 = (uint32) terminalData.out_stream.find('.', token.length());
                if (pos1 != std::string::npos) {
                    auto pos2 = (uint32) terminalData.out_stream.find('.', pos1+1);
                    // Ensure we were able to find last '.' and can use that index in string
                    if (pos2 != std::string::npos && pos2 + 1 < terminalData.out_stream.length()) {

                        // Successfully can parse version
                        mVersionMajor = (uint16) std::stoi(terminalData.out_stream.substr(token.length(), pos1 - token.length()));
                        mVersionMinor = (uint16) std::stoi(terminalData.out_stream.substr(pos1+1, pos2-1-pos1));
                        mVersionRev   = (uint16) std::stoi(terminalData.out_stream.substr(pos2+2, terminalData.out_stream.length()-pos2+1));

                        if (mVersionMajor == 4) {
                            mBuscoVersion = BUSCO_VERSION_4;
                        } else if (mVersionMajor == 3){
                            mBuscoVersion = BUSCO_VERSION_3;
                        } else {
                            FS_dprint("WARNING could not find version, assuming BUSCO version 4.0.2");
                            mBuscoVersion = BUSCO_VERSION_4;
                            ret = false;
                        }

                    } else {
                        FS_dprint("WARNING unable to find BUSCO version");
                        ret = false;
                    }
                } else {
                    FS_dprint("WARNING unable to find '.' in BUSCO version");
                    ret = false;
                }

            } else {
                FS_dprint("WARNING unable to parse BUSCO version");
                ret = false;
            }

        } else {
            FS_dprint("WARNING unable to execute BUSCO command to determine version");
            ret = false;
        }


    } catch (std::exception &e) {
        FS_dprint("BUSCO version standard exception: "+ std::string(e.what()));
        ret = false;
    } catch (...) {
        FS_dprint("BUSCO version unhandled exception...");
        ret = false;
    }

    if (!ret) {
        FS_dprint("WARNING unable to find BUSCO version, assuming defaults...");
        set_version_defaults();
    }

    FS_dprint("Success! BUSCO version set to: " + print_version());
    return ret;
}

void ModBUSCO::set_version_defaults() {
    FS_dprint("Setting BUSCO version defaults: " +
              std::to_string(DEFAULT_VERSION_MAJOR) + "." +
              std::to_string(DEFAULT_VERSION_MINOR) + "." +
              std::to_string(DEFAULT_VERSION_REV));
    mVersionMajor = DEFAULT_VERSION_MAJOR;
    mVersionMinor = DEFAULT_VERSION_MINOR;
    mVersionRev   = DEFAULT_VERSION_REV;
    mBuscoVersion = DEFAULT_VERSION;
}

bool ModBUSCO::is_version_supported(BUSCO_VERSION &version) {
    FS_dprint("Checking BUSCO version support...");

    bool ret = true;

    switch (version) {
        case BUSCO_VERSION_UNKNOWN:
            FS_dprint("WARNING Unknown version is NOT supported");
            ret = false;
            break;

        case BUSCO_VERSION_3:
            FS_dprint("WARNING BUSCO version 3 is NOT supported");
            ret = false;
            break;

        case BUSCO_VERSION_4:
            FS_dprint("BUSCO version 4 IS supported");
            break;
    }
    return ret;
}

std::string ModBUSCO::print_version() {
    std::string ret;

    ret =
          std::to_string(DEFAULT_VERSION_MAJOR) + "." +
          std::to_string(DEFAULT_VERSION_MINOR) + "." +
          std::to_string(DEFAULT_VERSION_REV);
    return ret;
}

std::string ModBUSCO::get_final_table_path(std::string &output_dir, std::string &database, BUSCO_VERSION &version) {
    // WARNING must have version by now

    std::string ret;

    switch (version) {

        case BUSCO_VERSION_UNKNOWN:
            break;

        case BUSCO_VERSION_3:
        case BUSCO_VERSION_4:
            // BUSCO Version 3 and 4 output the full_table.tsv as follows with eukaryota_odb1 database
            // OUTPUT_FLAG/run_eukaryota_odb10/full_table.tsv

            // This is the run_eukaryota_odb10 directory
            std::string busco_run_dir = BUSCO_PREPEND_TAG + mpFileSystem->get_filename(database, false);

            // Combine the run directory with the full_table.tsv file [run_eukaryota_odb10/full_table.tsv]
            std::string temp = PATHS(busco_run_dir, BUSCO_FULL_TABLE_FILENAME);

            // Combine final paths to current working directory
            // WARNING BUSCO 3,4 always output to CWD
            ret = PATHS(output_dir, temp);
            break;
    }

    return ret;
}
