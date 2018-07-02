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

ModEggnogDMND::ModEggnogDMND(std::string &out, std::string &in_hits, std::string &ont_out, bool blastp,
                             std::vector<uint16> &lvls, EntapDataPtrs &entap_data, std::string sql_db_path)
        : AbstractOntology(out, in_hits, ont_out, blastp, lvls, entap_data) {

    _egg_out_dir = PATHS(_ontology_dir, EGGNOG_DMND_DIR);
    _eggnog_db_path = sql_db_path;
}

std::pair<bool, std::string> ModEggnogDMND::verify_files() {
    FS_dprint("Overwrite was unselected, verifying output files...");
    _out_hits = get_output_dmnd_filepath();


    return std::pair<bool, std::string>();
}

void ModEggnogDMND::execute() {
    std::string                        std_out;
    std::string                        cmd;
    std::string                        blast;

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
    std_out = _out_hits + FileSystem::EXT_STD;

    //Run DIAMOND
    cmd =
            DIAMOND_EXE + " " +
            blast +
            " -d " + EGG_DMND_PATH +
            " --top 3"             +
            " --more-sensitive"    +
            " -q "                 + _inpath   +
            " -o "                 + _out_hits +
            " -p "                 + std::to_string(_threads) +
            " -f " + "6 qseqid sseqid pident length mismatch gapopen "
                    "qstart qend sstart send evalue bitscore qcovhsp stitle";

    if (TC_execute_cmd(cmd, std_out) != 0) {
        // Error in run
        _pFileSystem->delete_file(_out_hits);
        throw ExceptionHandler("Error in running DIAMOND against EggNOG database at: " +
                               EGG_DMND_PATH, ERR_ENTAP_RUN_EGGNOG);
    }
}

void ModEggnogDMND::parse() {

}

bool ModEggnogDMND::is_executable() {
//    std::string test_cmd;
//    uint8       err_code;

    return false;
}


std::string ModEggnogDMND::get_output_dmnd_filepath() {
    std::string filename;

    _blastp ? filename = "blastp" : filename = "blastx";
    filename += "_" + _pUserInput->get_user_transc_basename() + "_eggnog_db" + FileSystem::EXT_OUT;
    return PATHS(_egg_out_dir, filename);
}
