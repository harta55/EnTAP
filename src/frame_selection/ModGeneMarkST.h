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


#ifndef ENTAP_MODGENEMARKST_H
#define ENTAP_MODGENEMARKST_H


#include "AbstractFrame.h"

/**
 * ======================================================================
 * @class ModGeneMarkST
 *
 * Description          - This EnTAP module supports execution, parsing, and
 *                        statistical analysis of the GeneMarkS-T software
 *                        through terminal commands
 *                      - GeneMarkS-T performs Frame Selection on the
 *                        transcriptome input from the user. Sequences
 *                        where a frame was not found will be removed from
 *                        further pipeline analysis
 *                      - Parsed data is added to QueryData class
 *                      - Inherits from AbstractFrame and EntapModule classes
 *
 * Citation             - Shiyuyun Tang, Alexandre Lomsadze, Mark Borodovsky,
 *                        Identification of protein coding regions in RNA
 *                        transcripts, Nucleic Acids Research, Volume 43,
 *                        Issue 12, 13 July 2015, Page e78,
 *                        https://doi.org/10.1093/nar/gkv227
 *                      - http://exon.gatech.edu/GeneMark/
 *
 * ======================================================================
 */
class ModGeneMarkST : public AbstractFrame{

// Data pulled from GeneMarkS-T FASTA files
struct FrameData {
    std::string sequence;   // Sequence (either nucleotide protein) WITH sequence ID
    std::string frame_type; // Frame type (partial, removed, 5', 3', complete)
};
typedef std::map<std::string,ModGeneMarkST::FrameData> frame_map_t;

public:

    //******************* Public Functions *********************
    ModGeneMarkST(std::string &execution_stage_path, std::string &in_hits,
                  EntapDataPtrs &entap_data);
    ~ModGeneMarkST();
    virtual ModVerifyData verify_files() override ;
    virtual void execute() override ;
    virtual void parse() override ;
    virtual bool set_version() override;
    static bool is_executable(std::string& exe);
    //**********************************************************

private:

    typedef enum {

        GENEMARK_RETURN_OK=1        // WARNING perl exit commands need to be shifted by 8 >>8

    } GENEMARK_RETURN_CODE;

    //****************** Private Functions *********************
    void genemark_parse_fasta(std::string& fasta, FileSystem::ENT_FILE_TYPES file_type);
    void genemark_parse_lst(std::string & lst_path);
    //**********************************************************

    //**************** Private Const Variables *****************
    const std::string GENEMARK_LOG_FILE             = "gms.log";
    const std::string GENEMARK_HMM_FILE             = "GeneMark_hmm.mod";
    const std::string GENEMARK_STD_OUT              = "genemark_run";
    static std::vector<ENTAP_HEADERS> DEFAULT_HEADERS;

    //**********************************************************

    //****************** Private Variables *********************
    std::string mFinalFnnPath;      // Absolute path to FNN (nucleotide) file produced from GeneMarkS-T
    std::string mFinalLstPath;      // Absolute path to .lst file produced from GeneMarkS-T
    std::string mFinalGmstLogPath;  // Absolute path to gmst file produced from GeneMarkS-T
    std::string mFinalHmmPath;      // Absolute path to HMM file produced from GeneMarkS-T
    std::string mTranscriptomeFilename;    // Filename of input transcriptome
    //**********************************************************
};


#endif //ENTAP_MODGENEMARKST_H
