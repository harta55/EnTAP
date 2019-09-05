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

    mNewInput          = input;
    mpQueryData         = entap_data.mpQueryData;
    mpFileSystem        = entap_data.mpFileSystem;
    mpUserInput         = entap_data.mpUserInput;

    mEntapDataPtrs = entap_data;

    mOutpath            = mpFileSystem->get_root_path();
    mIsOverwrite       = mpUserInput->has_input(mpUserInput->INPUT_FLAG_OVERWRITE);
    mSoftwareFlags     = mpUserInput->get_user_input<vect_uint16_t>(mpUserInput->INPUT_FLAG_ONTOLOGY);
    mGoLevels          = mpUserInput->get_user_input<vect_uint16_t>(mpUserInput->INPUT_FLAG_GO_LEVELS);
    mOntologyDir       = PATHS(mOutpath, ONTOLOGY_OUT_PATH);
    mFinalOutputDir  = mpFileSystem->get_final_outdir();
    mEggnogDbPath     = EGG_SQL_DB_PATH;
    mInterproDatabases = mpUserInput->get_user_input<vect_str_t>(mpUserInput->INPUT_FLAG_INTERPRO);
    mAlignmentFileTypes = mpUserInput->get_user_output_types();

    if (mIsOverwrite) mpFileSystem->delete_dir(mOntologyDir);
    mpFileSystem->create_dir(mOntologyDir);
    mpFileSystem->delete_dir(mFinalOutputDir);
    mpFileSystem->create_dir(mFinalOutputDir);
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
        for (uint16 software : mSoftwareFlags) {
            ptr = spawn_object(software);
            verify_data = ptr->verify_files();
            if (!verify_data.files_exist) ptr->execute();
            ptr->parse();
            ptr.reset();
        }
        // If no error, flag
        mpQueryData->set_is_success_ontology(true);

        // Print final annotations
        print_eggnog();
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
                    mOutpath,
                    mNewInput,
                    mOntologyDir,
                    mBlastp,
                    mGoLevels,
                    mEntapDataPtrs,
                    mEggnogDbPath    // additional data
            ));
#endif
        case ONT_INTERPRO_SCAN:
            return std::unique_ptr<AbstractOntology>(new ModInterpro(
                    mOntologyDir,
                    mNewInput,
                    mEntapDataPtrs,
                    INTERPRO_EXE,
                    mInterproDatabases       // Additional data
            ));
        case ONT_EGGNOG_DMND:
            return std::unique_ptr<AbstractOntology>(new ModEggnogDMND(
                    mOntologyDir,
                    mNewInput,
                    mEntapDataPtrs,
                    DIAMOND_EXE,
                    mEggnogDbPath
            ));
        default:
            return std::unique_ptr<AbstractOntology>(new ModEggnogDMND(
                    mOntologyDir,
                    mNewInput,
                    mEntapDataPtrs,
                    DIAMOND_EXE,
                    mEggnogDbPath
            ));
    }
}


/**
 * ======================================================================
 * Function void Ontology::print_eggnog()
 *
 * Description          - Handles printing of final annotation output
 *                      - Current prints tsv file for all go levels specified,
 *                        no contam + contam files
 *
 * Notes                - None
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void Ontology::print_eggnog() {
    FS_dprint("Beginning to print final results...");

    std::string final_annotations_base;
    std::string final_annotations_contam_base;
    std::string final_annotations_no_contam_base;

    // TODO move to QueryData
    // Create output files for go levels (contaminants, no contam, all) and write headers
    for (uint16 lvl : mGoLevels) {

        final_annotations_base             = PATHS(mFinalOutputDir, FINAL_ANNOT_FILE);
        final_annotations_contam_base      = PATHS(mFinalOutputDir, FINAL_ANNOT_FILE_CONTAM);
        final_annotations_no_contam_base   = PATHS(mFinalOutputDir, FINAL_ANNOT_FILE_NO_CONTAM);

        mpQueryData->start_alignment_files(final_annotations_base, _HEADERS, (uint8)lvl, mAlignmentFileTypes);
        mpQueryData->start_alignment_files(final_annotations_contam_base, _HEADERS, (uint8)lvl, mAlignmentFileTypes);
        mpQueryData->start_alignment_files(final_annotations_no_contam_base, _HEADERS,(uint8) lvl, mAlignmentFileTypes);

        for (auto &pair : *mpQueryData->get_sequences_ptr()) {
            mpQueryData->add_alignment_data(final_annotations_base, pair.second, nullptr);

            if (pair.second->is_contaminant()) {
                mpQueryData->add_alignment_data(final_annotations_contam_base, pair.second, nullptr);
            } else {
                mpQueryData->add_alignment_data(final_annotations_no_contam_base, pair.second, nullptr);
            }
        }

        mpQueryData->end_alignment_files(final_annotations_base);
        mpQueryData->end_alignment_files(final_annotations_contam_base);
        mpQueryData->end_alignment_files(final_annotations_no_contam_base);
    }
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
    for (uint16 &flag : mSoftwareFlags) {
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