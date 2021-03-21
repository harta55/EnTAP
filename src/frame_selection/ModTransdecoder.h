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

#ifndef ENTAP_MODTRANSDECODER_H
#define ENTAP_MODTRANSDECODER_H


#include "AbstractFrame.h"

/**
 * ======================================================================
 * @class ModGeneMarkST
 *
 * Description          - This EnTAP module supports execution, parsing, and
 *                        statistical analysis of the TransDecoder software
 *                        through terminal commands
 *                      - TransDecoder performs Frame Selection by means of training
 *                        on the transcriptome input from the user. Sequences
 *                        where a frame was not found will be removed from
 *                        further pipeline analysis
 *                      - Parsed data is added to QueryData class
 *                      - Inherits from AbstractFrame and EntapModule classes
 *
 * Citation             - Haas, B. J., Papanicolaou, A., Yassour, M., Grabherr, M.,
 *                        Blood, P. D., Bowden, J., … Regev, A. (2013). De novo
 *                        transcript sequence reconstruction from RNA-seq using
 *                        the Trinity platform for reference generation and analysis.
 *                        Nature protocols, 8(8), 1494–1512. doi:10.1038/nprot.2013.084
 *                      - https://github.com/TransDecoder/TransDecoder/wiki
 *
 * ======================================================================
 */
class ModTransdecoder : public AbstractFrame {

public:

    //******************* Public Functions *********************
    ModTransdecoder(std::string &execution_stage_path, std::string &in_hits,
            EntapDataPtrs &entap_data);
    ~ModTransdecoder();

    virtual ModVerifyData verify_files() override;
    virtual void execute() override ;
    virtual void parse() override;
    virtual bool set_version() override;

    static bool is_executable(std::string &long_orfs_exe, std::string &predict_exe);
    //**********************************************************

private:

    //****************** Private Functions *********************
    int train_data(std::string &err_msg);
    int predict_frame(std::string &err_msg);
    void parse_transdecoder_fasta(std::string &file_path, FileSystem::ENT_FILE_TYPES file_type);
    void parse_transdecoder_fasta_header(std::string &seq_id, std::string &line, std::string &header, fp32 &score);
    bool verify_ORF(std::string& frame);
    //**********************************************************


    //**************** Private Const Variables *****************

    /* Transdecoder Commands as of v5.5.0 */
    const std::string CMD_TRANSCRIPTOME_INPUT = "-t";      // Command to specify input
    const std::string CMD_MIN_PROTEIN_LENGTH  = "-m";
    const std::string CMD_OUTPUT_DIR          = "-O";                   // Specify output directory WARNING added in v5.5.0
                                                                        // and doesn't really work
    const std::string CMD_SINGLE_BEST_ONLY    = "--single_best_only";   // Only retain one best ORF per transcript
    const std::string CMD_NO_REFINE_STARTS    = "--no_refine_starts";   // start refinement identifies potential start codons for 5' partial
                                                                        // ORFs using a PWM, process on by default.

    const std::string STD_OUTPUT_TRAINING = "long_orfs_std";
    const std::string STD_OUTPUT_PREDICTION = "prediction_std";

    // String flags in cds/pep output files to indicate FRAME
    const std::string FLAG_FRAME_BEGIN = "type:";                   // Marks beginning of frame type in output files
    const char        FLAG_FRAME_END   = ' ';
    const std::string FLAG_SCORE_BEGIN = "score=";
    const char        FLAG_SCORE_END   = ' ';
    const std::string TRANSDECODER_FLAG_COMPLETE_FRAME = "complete";
    const std::string TRANSDECODER_FLAG_INTERNAL_FRAME = "internal";
    const std::string TRANSDECODER_FLAG_5PRIME_FRAME = "5prime_partial";
    const std::string TRANSDECODER_FLAG_3PRIME_FRAME = "3prime_partial";
    // String flags in cds/pep output files surrounding QUERY ID
    // (ex: >MSTRG.1001.1.p1 GENE.MSTRG.1001.1~~MSTRG.1001.1.p1)
    const std::string FLAG_QUERYID_BEGIN = "GENE.";
    const std::string FLAG_QUERYID_END   = "~~";
    // Sometimes, TransDecoder can report the following:
    //  >TRINITY_DN0_c1_g1_i3.p1 TRINITY_DN0_c1_g1~~TRINITY_DN0_c1_g1_i3.p1
    //  ORF type:internal len:121 (-),score=46.75 TRINITY_DN0_c1_g1_i3:3-362(-)
    const char FLAG_QUERYID_V2_BEGIN = ' ';  // Space after 'score:'
    const char FLAG_QUERYID_V2_END = ':';   // some risk in this if the user has '.p' in sequence
    // Transdecoder file naming conventions
    const std::string FILE_TRANSDECODER_SUFFIX = ".transdecoder";

    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;
    //**********************************************************


    //****************** Private Variables *********************
    std::string mTransdecoderLongOrfsExe;       // Method of execution for LongOrds (train) portion of Transdecoder from user
    std::string mTransdecoderPredictExe;        // Method of execution for Prediction portion of Transdecoder from user
    std::string mOutputCDSFilePath;             // Absolute path to CDS file output from TransDecoder (nucleo)
    bool        mIsNoRefineStarts;              // User wants to pipe command '--no_refine_starts' to TransDecoder
    uint16 mMinProteinLength;
    //**********************************************************
};


#endif //ENTAP_MODTRANSDECODER_H
