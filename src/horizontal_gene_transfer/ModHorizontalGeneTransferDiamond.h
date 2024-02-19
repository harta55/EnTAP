/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2023, Alexander Hart, Dr. Jill Wegrzyn
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

#ifndef ENTAP_MODHORIZONTALGENETRANSFERDIAMOND_H
#define ENTAP_MODHORIZONTALGENETRANSFERDIAMOND_H


#include "AbstractHorizontalGeneTransfer.h"
#include "../common.h"
#include "../similarity_search/AbstractSimilaritySearch.h"
#include "../similarity_search/ModDiamond.h"
#include "../ExceptionHandler.h"
#include "../GraphingManager.h"

// TODO restructure this to be more modular with DIAMOND module and overhaul DIAMOND mod
//  there's a lot of duplicate code between this and that module that will need to be overhauled for
class ModHorizontalGeneTransferDiamond : public AbstractHorizontalGeneTransfer {

public:
    ModHorizontalGeneTransferDiamond(std::string &execution_stage_path, std::string &fasta_path, EntapDataPtrs &entap_data);
    ~ModHorizontalGeneTransferDiamond()=default;

    // ModEntap overrides
    virtual ModVerifyData verify_files() override;
    virtual void execute() override ;
    virtual void parse() override ;
    static bool is_executable(std::string& exe);
    virtual bool set_version() override;

    enum HGT_DATABASE_TYPES {
        HGT_DATABASE_UNKNOWN=0,
        HGT_DATABASE_DONOR,
        HGT_DATABASE_RECIPIENT
    };
protected:


    struct HGTDatabase {
        std::string database_path;      // Full path to database
        std::string database_shortname; // Shortname of the database (basically filename)
        std::string diamond_output;     // Output path after DIAMOND has been ran

        bool diamond_ran_success;       // Has DIAMOND been ran against this database?
        HGT_DATABASE_TYPES database_type;
    };
    static constexpr int DMND_COL_NUMBER = 14;

    // Terminal Commands (as of DIAMOND v0.9.9)
    // WARNING until restructuring of code make sure this matches ModDiamond.h
    const std::string CMD_QUERY_COVERAGE   = "--query-cover";     // Specify minimum query coverage for alignment
    const std::string CMD_SUBJECT_COVERAGE = "--subject-cover";   // Specify minimum target coverage for alignment
    const std::string CMD_MORE_SENSITIVE   = "--very-sensitive";  // Specify 'very sensitive' run that will take longer
    const std::string CMD_EVALUE           = "--evalue";          // Specify highest e-value to accept alignments for
    const std::string CMD_BLASTX           = "blastx";
    const std::string CMD_BLASTP           = "blastp";
    const std::string CMD_DATABASE         = "-d";                // Target database to align against
    const std::string CMD_QUERY_PATH       = "-q";                // Path to Query FASTA file
    const std::string CMD_OUTPUT_PATH      = "-o";                // Path to output
    const std::string CMD_THREADS          = "-p";                // Number of threads to use
    const std::string CMD_TOP_ALIGNMENTS   = "--max-target-seqs"; // Only keep top alignments (integer)
    const std::string CMD_OUTPUT_FORMAT    = "-f";
    //**********************************************************
    // Terminal Command EntapDefaults
    const uint16 CMD_DEFAULT_TOP_ALIGN  = 3;
    const std::string CMD_DEFAULT_OUTPUT_FORMAT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle";

    // Stats output defines
    const std::string DIAMOND_PREFIX                             = "hgt_diamond_";  // Prefix to be used when printing out files
    const std::string HGT_DIAMOND_DATABASE_BEST_HITS              = "annotated";
    const std::string HGT_DIAMOND_DATABASE_NO_HITS                = "unannotated";
    const std::string HGT_DIAMOND_DATABASE_UNSELECTED             = "unselected_hits";
    const std::string UNINFORMATIVE_FLAG                         = "Uninformative";
    const std::string INFORMATIVE_FLAG                           = "Informative";
    const std::string NO_HIT_FLAG                                = "No Hits";

    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;
    std::vector<HGTDatabase> mHGTDatabases;

    void calculate_best_stats (ModHorizontalGeneTransferDiamond::HGTDatabase &database);
    void calculate_hgt_candidates(std::vector<HGTDatabase> hgt_databases);
    EntapModule::ModVerifyData verify_databases(ent_input_multi_str_t &databases, HGT_DATABASE_TYPES data_type);
    bool run_diamond(AbstractSimilaritySearch::SimSearchCmd *cmd);

};


#endif //ENTAP_MODHORIZONTALGENETRANSFERDIAMOND_H
