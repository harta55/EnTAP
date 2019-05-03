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

#include "Ontology.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
#include "ontology/AbstractOntology.h"
#include "ontology/ModEggnog.h"
#include "ontology/ModInterpro.h"
#include "FileSystem.h"
#include "ontology/ModEggnogDMND.h"
#include "similarity_search/ModDiamond.h"

const int16 FINAL_ANNOT_LEN = 3;


/**
 * ======================================================================
 * Function Ontology::Ontology(int thread, std::string outpath, std::string input,
                   boost::program_options::variables_map &user_input, GraphingManager* graphing,
                   QueryData *queryData, bool blastp)
 *
 * Description          - Initializes Ontology object with values from user
 *                        input map/pushed vals
 *                      - Receives info on software that will be ran
 *
 * Notes                - Constructor
 *
 * @param thread        - Thread count
 * @param outpath       - Output directory specified
 * @param input         - Input fasta file to use (overwritten on execute if normal state)
 * @param user_input    - Boost map of user input
 * @param graphing      - Ptr to graphing manager instance
 * @param queryData     - Ptr to queryData instance
 * @param blastp        - Blastp yes/no
 *
 * @return              - None
 *
 * =====================================================================
 */
Ontology::Ontology(std::string input, EntapDataPtrs &entap_data) {
    FS_dprint("Spawn Object - Ontology");

    _new_input          = input;
    _pGraphingManager   = entap_data._pGraphingManager;
    _pQueryData         = entap_data._pQueryData;
    _pFileSystem        = entap_data._pFileSystem;
    _pUserInput         = entap_data._pUserInput;
    _pEntapDatabase     = entap_data._pEntapDatbase;

    _entap_data_ptrs = entap_data;

    _threads            = _pUserInput->get_supported_threads();
    _outpath            = _pFileSystem->get_root_path();
    _is_overwrite       = _pUserInput->has_input(_pUserInput->INPUT_FLAG_OVERWRITE);
    _software_flags     = _pUserInput->get_user_input<vect_uint16_t>(_pUserInput->INPUT_FLAG_ONTOLOGY);
    _go_levels          = _pUserInput->get_user_input<vect_uint16_t>(_pUserInput->INPUT_FLAG_GO_LEVELS);
    _blastp             = _pUserInput->has_input(_pUserInput->INPUT_FLAG_RUNPROTEIN);
    _ontology_dir       = PATHS(_outpath, ONTOLOGY_OUT_PATH);
    _final_outpath_dir  = _pFileSystem->get_final_outdir();
    _eggnog_db_path     = EGG_SQL_DB_PATH;
    _interpro_databases = _pUserInput->get_user_input<vect_str_t>(_pUserInput->INPUT_FLAG_INTERPRO);
    _alignment_file_types = _pUserInput->get_user_output_types();

    if (_is_overwrite) _pFileSystem->delete_dir(_ontology_dir);
    _pFileSystem->create_dir(_ontology_dir);
    _pFileSystem->delete_dir(_final_outpath_dir);
    _pFileSystem->create_dir(_final_outpath_dir);
}


/**
 * ======================================================================
 * Function void Ontology::execute(std::string input,std::string no_hit)
 *
 * Description          - Manager of running/parsing software to be ran
 *                        for ontology analysis
 *                      - Runs separate analyses for hits and no hits from sim search
 *
 * Notes                - Execution entry
 *
 * @param input         - Fasta sequences that hit databases
 * @param no_hit        - Fasta sequences that did not hit sequences
 *
 * @return              - None
 *
 * =====================================================================
 */
void Ontology::execute() {

    EntapModule::ModVerifyData verify_data;
    std::unique_ptr<EntapModule> ptr;

    init_headers();
    try {
        for (uint16 software : _software_flags) {
            ptr = spawn_object(software);
            verify_data = ptr->verify_files();
            if (!verify_data.files_exist) ptr->execute();
            ptr->parse();
            ptr.reset();
        }
        print_eggnog(*_pQueryData->get_sequences_ptr());
    } catch (ExceptionHandler &e) {
        ptr.reset();
        throw e;
    }
}


/**
 * ======================================================================
 * Function std::unique_ptr<AbstractOntology> Ontology::spawn_object(uint16 &software)
 *
 * Description          - Spawns object for specified ontology software
 *
 * Notes                - None
 *
 * @param software      - Int flag to specify software
 *
 * @return              - Ptr to ontology software object
 *
 * =====================================================================
 */
std::unique_ptr<AbstractOntology> Ontology::spawn_object(uint16 &software) {
    switch (software) {
#ifdef EGGNOG_MAPPER
        case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
            return std::unique_ptr<AbstractOntology>(new ModEggnog(
                    _outpath,
                    _new_input,
                    _ontology_dir,
                    _blastp,
                    _go_levels,
                    _entap_data_ptrs,
                    _eggnog_db_path    // additional data
            ));
#endif
        case ONT_INTERPRO_SCAN:
            return std::unique_ptr<AbstractOntology>(new ModInterpro(
                    _ontology_dir,
                    _new_input,
                    _entap_data_ptrs,
                    INTERPRO_EXE,
                    _interpro_databases       // Additional data
            ));
        case ONT_EGGNOG_DMND:
            return std::unique_ptr<AbstractOntology>(new ModEggnogDMND(
                    _ontology_dir,
                    _new_input,
                    _entap_data_ptrs,
                    DIAMOND_EXE,
                    _eggnog_db_path
            ));
        default:
            return std::unique_ptr<AbstractOntology>(new ModEggnogDMND(
                    _ontology_dir,
                    _new_input,
                    _entap_data_ptrs,
                    DIAMOND_EXE,
                    _eggnog_db_path
            ));
    }
}


/**
 * ======================================================================
 * Function void Ontology::print_eggnog(QUERY_MAP_T &SEQUENCES)
 *
 * Description          - Handles printing of final annotation output
 *                      - Current prints tsv file for all go levels specified,
 *                        no contam + contam files
 *
 * Notes                - None
 *
 * @param SEQUENCES     - Map of sequence data
 *
 * @return              - None
 *
 * =====================================================================
 */
void Ontology::print_eggnog(QUERY_MAP_T &SEQUENCES) {
    FS_dprint("Beginning to print final results...");
    std::map<uint16, std::ofstream*[FINAL_ANNOT_LEN]> file_map;
    std::string file_name;
    std::string file_contam;
    std::string file_no_contam;
    std::string outpath;
    std::string out_contam;
    std::string out_no_contam;

    std::string final_annotations_base;
    std::string final_annotations_contam_base;
    std::string final_annotations_no_contam_base;

    // TODO move to QueryData
    // Create output files for go levels (contaminants, no contam, all) and write headers
    for (uint16 lvl : _go_levels) {

        final_annotations_base             = PATHS(_final_outpath_dir, FINAL_ANNOT_FILE);
        final_annotations_contam_base      = PATHS(_final_outpath_dir, FINAL_ANNOT_FILE_CONTAM);
        final_annotations_no_contam_base   = PATHS(_final_outpath_dir, FINAL_ANNOT_FILE_NO_CONTAM);

        _pQueryData->start_alignment_files(final_annotations_base, _HEADERS, (uint8)lvl, _alignment_file_types);
        _pQueryData->start_alignment_files(final_annotations_contam_base, _HEADERS, (uint8)lvl, _alignment_file_types);
        _pQueryData->start_alignment_files(final_annotations_no_contam_base, _HEADERS,(uint8) lvl, _alignment_file_types);

        for (auto &pair : SEQUENCES) {
            _pQueryData->add_alignment_data(final_annotations_base, pair.second);

            if (pair.second->isContaminant()) {
                _pQueryData->add_alignment_data(final_annotations_contam_base, pair.second);
            } else {
                _pQueryData->add_alignment_data(final_annotations_no_contam_base, pair.second);
            }
        }

        _pQueryData->end_alignment_files(final_annotations_base);
        _pQueryData->end_alignment_files(final_annotations_contam_base);
        _pQueryData->end_alignment_files(final_annotations_no_contam_base);
    }
//
//    // Write to all files
//    for (auto &pair : SEQUENCES) {
//        for (uint16 lvl : _go_levels) {
//
//            for (uint16 i=0; i < FINAL_ANNOT_LEN; i++) {
//                if (i == FINAL_ALL_IND) {
//                    *file_map[lvl][i] << pair.second->print_delim(_HEADERS, lvl, FileSystem::DELIM_TSV) << std::endl;
//                } else if (i == FINAL_CONTAM_IND && pair.second->isContaminant()) {
//                    *file_map[lvl][i]<< pair.second->print_delim(_HEADERS, lvl, FileSystem::DELIM_TSV)<<std::endl;
//                } else if (i == FINAL_NO_CONTAM_IND && !pair.second->isContaminant()) {
//                    *file_map[lvl][i]<< pair.second->print_delim(_HEADERS, lvl, FileSystem::DELIM_TSV)<<std::endl;
//                }
//            }
//        }
//    }
//
//    // Cleanup
//    for(auto& pair : file_map) {
//        for (uint16 i=0; i < FINAL_ANNOT_LEN; i++) {
//            pair.second[i]->close();
//            delete pair.second[i];
//        }
//    }
    FS_dprint("Success!");
}


/**
 * ======================================================================
 * Function void Ontology::init_headers()
 *
 * Description          - Initializes default headers that will be printed
 *                        to final tsv as well as extra headers for software
 *                        being ran
 *
 * Notes                - None
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void Ontology::init_headers() {

    std::vector<ENTAP_HEADERS>     out_header;
    std::vector<ENTAP_HEADERS>     add_header;
    // Add default sim search headers (pulled from SimilaritySearch.c, separate in case we want something else)
    out_header = ModDiamond::DEFAULT_HEADERS;

    // Add additional headers for ontology software
    for (uint16 &flag : _software_flags) {
        switch (flag) {
#ifdef EGGNOG_MAPPER
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
#endif
            case ONT_INTERPRO_SCAN:
                add_header = {
                        ModInterpro::DEFAULT_HEADERS
                };
                break;
            case ONT_EGGNOG_DMND:
                add_header = {
                        ModEggnogDMND::DEFAULT_HEADERS
                };
                break;
            default:
                throw ExceptionHandler("ERROR: Unknown INT flag used during Ontology",
                    ERR_ENTAP_MEM_ALLOC);
        }
        out_header.insert(out_header.end(), add_header.begin(), add_header.end());
    }

    _HEADERS = out_header;
}