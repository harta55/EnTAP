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

#include "Ontology.h"
#include "ExceptionHandler.h"
#include "EntapGlobals.h"
#include "ontology/AbstractOntology.h"
#include "ontology/ModEggnog.h"
#include "ontology/ModInterpro.h"
#include "FileSystem.h"
#include "ontology/ModEggnogDMND.h"
#include "similarity_search/ModDiamond.h"
#include "ontology/ModBUSCO.h"

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

    mNewInput           = input;
    mpQueryData         = entap_data.mpQueryData;
    mpFileSystem        = entap_data.mpFileSystem;
    mpUserInput         = entap_data.mpUserInput;

    mEntapDataPtrs = entap_data;

    mOutpath            = mpFileSystem->get_root_path();
    mIsOverwrite        = mpUserInput->has_input(INPUT_FLAG_OVERWRITE);
    mSoftwareFlags      = mpUserInput->get_user_input<ent_input_multi_int_t>(INPUT_FLAG_ONTOLOGY);
    mGoLevels           = mpUserInput->get_user_input<ent_input_multi_int_t>(INPUT_FLAG_GO_LEVELS);
    mOntologyDir        = PATHS(mOutpath, ONTOLOGY_OUT_PATH);
    mFinalOutputDir     = mpFileSystem->get_final_outdir();
    mAlignmentFileTypes = mpUserInput->get_user_output_types();
    mEntapHeaders       = mpUserInput->get_user_input<std::vector<ENTAP_HEADERS>>(INPUT_FLAG_ENTAP_HEADERS);

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

    try {
        for (uint16 software : mSoftwareFlags) {
            ptr = spawn_object(software);
            verify_data = ptr->verify_files();
            if (!verify_data.files_exist) ptr->execute();
            ptr->parse();
            ptr->set_success_flags();
            ptr.reset();
        }

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
    ent_input_str_t exe_path;
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
                    mEntapDataPtrs
            ));
        case ONT_EGGNOG_DMND:
            return std::unique_ptr<AbstractOntology>(new ModEggnogDMND(
                    mOntologyDir,
                    mNewInput,
                    mEntapDataPtrs
            ));
        case ONT_BUSCO:
            return std::unique_ptr<AbstractOntology>(new ModBUSCO(
                    mOntologyDir,
                    mNewInput,
                    mEntapDataPtrs
             ));
        default:
            return std::unique_ptr<AbstractOntology>(new ModEggnogDMND(
                    mOntologyDir,
                    mNewInput,
                    mEntapDataPtrs
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

    final_annotations_base             = PATHS(mFinalOutputDir, FINAL_ANNOT_FILE);
    mpQueryData->start_alignment_files(final_annotations_base, mEntapHeaders, mGoLevels, mAlignmentFileTypes);

    final_annotations_contam_base      = PATHS(mFinalOutputDir, FINAL_ANNOT_FILE_CONTAM);
    mpQueryData->start_alignment_files(final_annotations_contam_base, mEntapHeaders, mGoLevels, mAlignmentFileTypes);

    final_annotations_no_contam_base   = PATHS(mFinalOutputDir, FINAL_ANNOT_FILE_NO_CONTAM);
    mpQueryData->start_alignment_files(final_annotations_no_contam_base, mEntapHeaders, mGoLevels, mAlignmentFileTypes);

    // Add respective alignements
    for (auto &pair : *mpQueryData->get_sequences_ptr()) {
        mpQueryData->add_alignment_data(final_annotations_base, pair.second, nullptr);

        if (pair.second->is_contaminant()) {
            mpQueryData->add_alignment_data(final_annotations_contam_base, pair.second, nullptr);
        } else {
            mpQueryData->add_alignment_data(final_annotations_no_contam_base, pair.second, nullptr);
        }
    }

    // Cleanup
    mpQueryData->end_alignment_files(final_annotations_base);
    mpQueryData->end_alignment_files(final_annotations_contam_base);
    mpQueryData->end_alignment_files(final_annotations_no_contam_base);
    FS_dprint("Success!");
}