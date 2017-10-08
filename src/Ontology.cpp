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

#include <boost/filesystem.hpp>
#include <csv.h>
#include <fstream>
#include <boost/archive/binary_iarchive.hpp>
#include "Ontology.h"
#include "EntapConfig.h"
#include <boost/serialization/map.hpp>
#include <iomanip>
#include "ExceptionHandler.h"
#include "EntapExecute.h"
#include "SimilaritySearch.h"
#include "DatabaseHelper.h"
#include "EntapGlobals.h"
#include "ontology/AbstractOntology.h"
#include "ontology/ModEggnog.h"
#include "ontology/ModInterpro.h"

Ontology::Ontology(int thread, std::string outpath, std::string input,
                   boost::program_options::variables_map &user_input, GraphingManager* graphing,
                   QueryData *queryData, bool blastp) {
    print_debug("Spawn object - Ontology");
    _ontology_exe       = EGG_EMAPPER_EXE;
    _threads            = (uint8)thread;
    _outpath            = outpath;
    _new_input          = input;
    _is_overwrite       = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _blastp             = blastp;
    _software_flags     = user_input[ENTAP_CONFIG::INPUT_FLAG_ONTOLOGY].as<std::vector<uint16>>();
    _go_levels          = user_input[ENTAP_CONFIG::INPUT_FLAG_GO_LEVELS].as<std::vector<uint16>>();
    _ontology_dir       = PATHS(outpath, ONTOLOGY_OUT_PATH);
    _eggnog_db_path     = EGG_SQL_DB_PATH;
    _graphingManager    = graphing;
    _QUERY_DATA         = queryData;
    _interpro_databases = user_input[ENTAP_CONFIG::INPUT_FLAG_INTERPRO].as<std::vector<std::string>>();
}


void Ontology::execute(std::string input,std::string no_hit) {

    std::pair<bool,std::string> verify_pair;
    std::unique_ptr<AbstractOntology> ptr;

    _new_input = input;
    _input_no_hits = no_hit;

    if (_is_overwrite) boostFS::remove_all(_ontology_dir);
    boostFS::create_directories(_ontology_dir);
    init_headers();
    try {
        for (uint16 software : _software_flags) {
            ptr = spawn_object(software);
            ptr->set_data(_eggnog_db_path, _interpro_databases);
            verify_pair = ptr->verify_files();
            if (!verify_pair.first) ptr->execute();
            ptr->parse();
            ptr.release();
        }
        print_eggnog(*_QUERY_DATA->get_sequences_ptr());
    } catch (ExceptionHandler &e) {throw e;}

    // TODO move printing to querydata
}


std::unique_ptr<AbstractOntology> Ontology::spawn_object(uint16 &software) {
    switch (software) {
        case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
            return std::unique_ptr<AbstractOntology>(new ModEggnog(
                    _ontology_exe, _outpath, _new_input, _input_no_hits,
                    _ontology_dir, _graphingManager, _QUERY_DATA, _blastp,
                    _go_levels, _threads
            ));
        case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
            break;
        default:
            return std::unique_ptr<AbstractOntology>(new ModInterpro(
                    _ontology_exe, _outpath, _new_input, _input_no_hits,
                    _ontology_dir, _graphingManager, _QUERY_DATA, _blastp,
                    _go_levels, _threads
            ));
    }
}


void Ontology::print_eggnog(query_map_struct &SEQUENCES) {
    print_debug("Beginning to print final results...");
    std::map<short, std::ofstream*> file_map;
    std::string file_name;
    std::string outpath;
    for (uint16 lvl : _go_levels) {
        file_name = "final_annotations_lvl" + std::to_string(lvl) + ".tsv";
        outpath = (boostFS::path(_outpath) / boostFS::path(file_name)).string();
        boostFS::remove(outpath);
        file_map[lvl] =
                new std::ofstream(outpath, std::ios::out | std::ios::app);
        for (const std::string *header : _HEADERS) {
            *file_map[lvl] << *header << '\t';
        }
        *file_map[lvl] << std::endl;
    }
    for (auto &pair : SEQUENCES) {
        for (uint16 i : _go_levels) {
            *file_map[i]<< pair.second.print_tsv(_HEADERS,i)<<std::endl;
        }
    }
    for(auto& pair : file_map) {
        pair.second->close();
        delete pair.second;
    }
    print_debug("Success!");
}


void Ontology::init_headers() {

    std::vector<const std::string*>     out_header;
    std::vector<const std::string*>     add_header;
    // Add default sim search headers
    out_header = {
            &ENTAP_EXECUTE::HEADER_QUERY,
            &ENTAP_EXECUTE::HEADER_SUBJECT,
            &ENTAP_EXECUTE::HEADER_PERCENT,
            &ENTAP_EXECUTE::HEADER_ALIGN_LEN,
            &ENTAP_EXECUTE::HEADER_MISMATCH,
            &ENTAP_EXECUTE::HEADER_GAP_OPEN,
            &ENTAP_EXECUTE::HEADER_QUERY_S,
            &ENTAP_EXECUTE::HEADER_QUERY_E,
            &ENTAP_EXECUTE::HEADER_SUBJ_S,
            &ENTAP_EXECUTE::HEADER_SUBJ_E,
            &ENTAP_EXECUTE::HEADER_E_VAL,
            &ENTAP_EXECUTE::HEADER_COVERAGE,
            &ENTAP_EXECUTE::HEADER_TITLE,
            &ENTAP_EXECUTE::HEADER_SPECIES,
            &ENTAP_EXECUTE::HEADER_DATABASE,
            &ENTAP_EXECUTE::HEADER_FRAME,
            &ENTAP_EXECUTE::HEADER_CONTAM,
            &ENTAP_EXECUTE::HEADER_INFORM
    };
    // Add additional headers for ontology software
    for (uint16 &flag : _software_flags) {
        switch (flag) {
            case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
                add_header = {
                        &ENTAP_EXECUTE::HEADER_SEED_ORTH,
                        &ENTAP_EXECUTE::HEADER_SEED_EVAL,
                        &ENTAP_EXECUTE::HEADER_SEED_SCORE,
                        &ENTAP_EXECUTE::HEADER_PRED_GENE,
                        &ENTAP_EXECUTE::HEADER_TAX_SCOPE,
                        &ENTAP_EXECUTE::HEADER_EGG_OGS,
                        &ENTAP_EXECUTE::HEADER_EGG_DESC,
                        &ENTAP_EXECUTE::HEADER_EGG_KEGG,
                        &ENTAP_EXECUTE::HEADER_EGG_PROTEIN,
                        &ENTAP_EXECUTE::HEADER_EGG_GO_BIO,
                        &ENTAP_EXECUTE::HEADER_EGG_GO_CELL,
                        &ENTAP_EXECUTE::HEADER_EGG_GO_MOLE
                };
                break;
            case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
                add_header = {
                        &ENTAP_EXECUTE::HEADER_INTER_EVAL,
                        &ENTAP_EXECUTE::HEADER_INTER_INTERPRO,
                        &ENTAP_EXECUTE::HEADER_INTER_DATA_TYPE,
                        &ENTAP_EXECUTE::HEADER_INTER_DATA_TERM,
                        &ENTAP_EXECUTE::HEADER_INTER_PATHWAY,
                        &ENTAP_EXECUTE::HEADER_INTER_GO_BIO,
                        &ENTAP_EXECUTE::HEADER_INTER_GO_CELL,
                        &ENTAP_EXECUTE::HEADER_INTER_GO_MOLE
                };
                break;
            default:
                break;
        }
        out_header.insert(out_header.end(), add_header.begin(), add_header.end());
    }
    _HEADERS = out_header;
}

void Ontology::print_header(std::string file) {
    std::ofstream ofstream(file, std::ios::out | std::ios::app);
    for (const std::string *val : _HEADERS) ofstream << *val << '\t';
    ofstream<<std::endl;
    ofstream.close();
}