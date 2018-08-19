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

#include "ModEggnogDMND.h"
#include "../readers/DelimReader.h"
#include "../database/EggnogDatabase.h"

ModEggnogDMND::ModEggnogDMND(std::string &out, std::string &in_hits, std::string &ont_out, bool blastp,
                             std::vector<uint16> &lvls, EntapDataPtrs &entap_data, std::string sql_db_path)
        : AbstractOntology(out, in_hits, ont_out, blastp, lvls, entap_data) {
    FS_dprint("Spawn Object - ModEggnogDMND");

    _egg_out_dir = PATHS(_ontology_dir, EGGNOG_DMND_DIR);
    _eggnog_db_path = sql_db_path;
}

std::pair<bool, std::string> ModEggnogDMND::verify_files() {
    bool file_exists = false;
    uint16 file_status = 0;

    FS_dprint("Overwrite was unselected, verifying output files...");
    _out_hits = get_output_dmnd_filepath(true);
    file_status = _pFileSystem->get_file_status(_out_hits);

    if (file_status != 0) {
        FS_dprint(_pFileSystem->print_file_status(file_status,_out_hits));
        FS_dprint("Errors in opening file, continuing with execution...");
    } else {
        file_exists = true;
    }
    return std::make_pair(file_exists, _out_hits);
}

void ModEggnogDMND::execute() {
    std::string                        std_out;
    std::string                        cmd;
    std::string                        blast;
    std::stringstream                  err_stream;
    std::stringstream                  out_stream;

    FS_dprint("Running EggNOG against Diamond database...");

    // Ensure both input path and EggNOG DMND database exist before continuing
    if (!_pFileSystem->file_exists(EGG_DMND_PATH)) {
        throw ExceptionHandler("EggNOG DIAMOND database not found at: " + EGG_DMND_PATH,
                               ERR_ENTAP_EGGNOG_FILES);
    }
    if (!_pFileSystem->file_exists(_inpath)) {
        throw ExceptionHandler("Input transcriptome not found at: " + _inpath, ERR_ENTAP_EGGNOG_FILES);
    }

    // Generate paths for DIAMOND run (out_hits set previously)
    std_out = get_output_dmnd_filepath(false) + "_" + FileSystem::EXT_STD;
    _blastp ? blast = "blastp" : blast = "blastx";

    //Run DIAMOND
    cmd =
            DIAMOND_EXE + " " +
            blast +
            " -d " + EGG_DMND_PATH +
            " --top 1"             +
            " --more-sensitive"    +
            " -q "                 + _inpath   +
            " -o "                 + _out_hits  +
            " -p "                 + std::to_string(_threads) +
            " -f " + "6 qseqid sseqid pident length mismatch gapopen "
                    "qstart qend sstart send evalue bitscore qcovhsp stitle";

    if (TC_execute_cmd(cmd, err_stream, out_stream, std_out) != 0) {
        // Error in run
        _pFileSystem->delete_file(_out_hits);
        FS_dprint("DIAMOND STD OUT:\n" + out_stream.str());
        throw ExceptionHandler("Error in running DIAMOND against EggNOG database at: " +
                               EGG_DMND_PATH + "\nDIAMOND Error:\n" + err_stream.str(), ERR_ENTAP_RUN_EGGNOG);
    }
}

void ModEggnogDMND::parse() {
    uint16 file_status=0;
    EggnogDatabase *eggnogDatabase;

    FS_dprint("Parsing EggNOG DMND file located at: " + _out_hits);

    // Ensure file is valid
    file_status = _pFileSystem->get_file_status(_out_hits);
    if (file_status != 0) {
        throw ExceptionHandler(_pFileSystem->print_file_status(file_status,_out_hits),
                               ERR_ENTAP_PARSE_EGGNOG_DMND);
    }

    // Generate EggNOG database
    eggnogDatabase = new EggnogDatabase(_pFileSystem);
    if (eggnogDatabase->open_sql(EGG_SQL_DB_PATH) != EggnogDatabase::ERR_EGG_OK) {
        throw ExceptionHandler("Unable to open EggNOG SQL Database", ERR_ENTAP_PARSE_EGGNOG_DMND);
    }

    // File valid, continue
#ifdef USE_FAST_CSV
    // ------------------ Read from DIAMOND output ---------------------- //
    std::string qseqid;
    std::string sseqid, stitle, database_name,pident, bitscore,
            length, mismatch, gapopen, qstart, qend, sstart, send;
    fp64 evalue;
    fp64 coverage;
    QuerySequence::EggnogResults eggnogResults;
    QuerySequence *querySequence;
    // ----------------------------------------------------------------- //
    // Begin using CSVReader lib to parse data
    io::CSVReader<DMND_COL_NUMBER, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(_out_hits);
    while (in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                       qstart, qend, sstart, send, evalue, bitscore, coverage,stitle)) {

        // Currently throwing away most DIAMOND results
        // Get SQL data
        eggnogResults = eggnogDatabase->get_eggnog_entry(sseqid);
        if (eggnogResults.member_ogs.empty()) {
            // Couldn't find a match
            FS_dprint("Couldn't find OGs for " + qseqid + " and " + sseqid);
        } else {
            // Found in database
            querySequence = _pQUERY_DATA->get_sequence(qseqid);
            if (querySequence == nullptr) {
                throw ExceptionHandler("Unable to find sequence " + qseqid + " in input transcriptome",
                    ERR_ENTAP_PARSE_EGGNOG_DMND);
            }
            querySequence->set_eggnog_results(eggnogResults);
        }
    }
#endif
}

bool ModEggnogDMND::is_executable() {
//    std::string test_cmd;
//    uint8       err_code;

    return false;
}

ModEggnogDMND::~ModEggnogDMND() {
    FS_dprint("Killing Object - ModEggnogDMND");

}


std::string ModEggnogDMND::get_output_dmnd_filepath(bool final) {
    std::string filename;

    _blastp ? filename = "blastp" : filename = "blastx";
    filename += "_" + _pUserInput->get_user_transc_basename() + "_eggnog_db";
    if (final) filename += FileSystem::EXT_OUT;
    return PATHS(_egg_out_dir, filename);
}
