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

#ifndef ENTAP_MODDIAMOND_H
#define ENTAP_MODDIAMOND_H


#include "AbstractSimilaritySearch.h"

/**
 * ======================================================================
 * @class ModInterpro
 *
 * Description          - This EnTAP module supports execution, parsing, and
 *                        statistical analysis of the DIAMOND software
 *                        through terminal commands
 *                      - DIAMOND will perform similarity searching against
 *                        DIAMOND compatible databases provided by the user
 *                      - Parsed data is added to QueryData class
 *                      - Inherits from AbstractSimilaritySearch and
 *                        EntapModule classes
 *
 * Citation             - P. Jones et al., “InterProScan 5: genome-scale
 *                        protein function classification,” (in eng),
 *                        Bioinformatics, vol. 30, no. 9, pp. 1236-40, May 2014.
 *
 * ======================================================================
 */
class ModDiamond : public AbstractSimilaritySearch {

public:
    //******************* Public Functions *********************
    ModDiamond(std::string &out, std::string &fasta_path,EntapDataPtrs &entap_data, vect_str_t &databases);
    ~ModDiamond() = default;

    // ModEntap overrides
    virtual ModVerifyData verify_files() override;
    virtual void execute() override ;
    virtual void parse() override ;
    static bool is_executable(std::string& exe);
    virtual bool set_version() override;

    // AbstractSimilaritySearch overrides
    bool run_blast(SimSearchCmd *cmd, bool use_defaults) override;
    //**********************************************************

private:
    //****************** Private Functions *********************
    void calculate_best_stats(bool is_final, std::string database_path="");
    //**********************************************************

    //************** Private Member Variables ******************
    bool mParsedFile;   // TRUE if any alignment file had alignments and was parsed
    //**********************************************************

    //**************** Private Const Variables *****************
    static constexpr int DMND_COL_NUMBER = 14;
    static std::vector<ENTAP_HEADERS> UNIPROT_HEADERS;
    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;

    // Graphing constants
    const uint8 GRAPH_SUM_FLAG                                   = 2;
    const std::string GRAPH_DATABASE_SUM_TITLE                   = "_Summary";
    const std::string GRAPH_DATABASE_SUM_TXT                     = "_summary_bar.txt";
    const std::string GRAPH_DATABASE_SUM_PNG                     = "_summary_bar.png";
    const std::string GRAPH_SPECIES_BAR_TXT                      = "_species_bar.txt";
    const std::string GRAPH_SPECIES_BAR_PNG                      = "_species_bar.png";
    const std::string GRAPH_SPECIES_TITLE                        = "_Top_10_Species_Distribution";
    const std::string GRAPH_CONTAM_BAR_TXT                       = "_contam_bar.txt";
    const std::string GRAPH_CONTAM_BAR_PNG                       = "_contam_bar.png";
    const std::string GRAPH_CONTAM_TITLE                         = "_Top_10_Contaminant_Distribution";
    const std::string UNINFORMATIVE_FLAG                         = "Uninformative";
    const std::string INFORMATIVE_FLAG                           = "Informative";
    const std::string NO_HIT_FLAG                                = "No Hits";

    // Terminal Command EntapDefaults
    const uint16 CMD_DEFAULT_TOP_ALIGN  = 3;
    const std::string CMD_DEFAULT_OUTPUT_FORMAT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle";

    // Terminal Commands (as of DIAMOND v0.9.9)
    const std::string CMD_QUERY_COVERAGE   = "--query-cover";     // Specify minimum query coverage for alignment
    const std::string CMD_SUBJECT_COVERAGE = "--subject-cover";   // Specify minimum target coverage for alignment
    const std::string CMD_MORE_SENSITIVE   = "--more-sensitive";  // Specify 'more sensitive' run that will take longer
    const std::string CMD_EVALUE           = "--evalue";          // Specify highest e-value to accept alignments for
    const std::string CMD_BLASTX           = "blastx";
    const std::string CMD_BLASTP           = "blastp";
    const std::string CMD_DATABASE         = "-d";                // Target database to align against
    const std::string CMD_QUERY_PATH       = "-q";                // Path to Query FASTA file
    const std::string CMD_OUTPUT_PATH      = "-o";                // Path to output
    const std::string CMD_THREADS          = "-p";                // Number of threads to use
    const std::string CMD_TOP_ALIGNMENTS   = "--top";             // Only keep top alignments (integer)
    const std::string CMD_OUTPUT_FORMAT    = "-f";
    //**********************************************************

    void set_uniprot_headers();
};


#endif //ENTAP_MODDIAMOND_H
