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
 * MERCHANTABILIT * but WITHOUT ANY WARRANTY; without even the implied warranty of
Y or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/


//*********************** Includes *****************************
#include "../common.h"
#include "ModGeneMarkST.h"
#include "../ExceptionHandler.h"
#include "../FileSystem.h"
#include "../TerminalCommands.h"
#include "../GraphingManager.h"
//**************************************************************


/**
 * ======================================================================
 * Function std::pair<bool, std::string> ModGeneMarkST::verify_files()
 *
 * Description           - Checks whether GeneMarkS-T has already been ran
 *                         with the same input
 *
 * Notes                 - None
 *
 *
 * @return               - Bool   - Yes if previous files have been found
 *                       - String - Path to frame selected transcriptome
 *
 * =====================================================================
 */
EntapModule::ModVerifyData ModGeneMarkST::verify_files() {
    ModVerifyData modVerifyData;    // Data to return from verification process

    modVerifyData.files_exist = false;

    FS_dprint("Beginning to verify GeneMarkS-T module files...");

    // If GeneMarkS-T has already been ran
    if (mpFileSystem->file_exists(mFinalFaaPath) && mpFileSystem->file_exists(mFinalLstPath)) {
        // Yes, files have been found and GeneMarkS-T has already been ran
        FS_dprint("Files found at: " + mFinalFaaPath + "\nand: " + mFinalLstPath +
                "\ncontinuing EnTAP with these files and skipping Frame Selection");
        modVerifyData.files_exist = true;
    } else {
        // No, files not found. GeneMarkS-T needs to be ran
        FS_dprint("File not found at " + mFinalFaaPath + "\nor " + mFinalLstPath +
                  " so continuing with Frame Selection");
    }
    modVerifyData.output_paths = vect_str_t{mFinalFaaPath};
    return modVerifyData;
}


/**
 * ======================================================================
 * Function std::string ModGeneMarkST::execute()
 *
 * Description          - Will run GeneMarkS-T through PStreams terminal
 *                        commands
 *                      - 4/5 output files are expected (3 are required to continue to parsing)
 *
 * Notes                - Fatal Error if an issue occurred during GeneMarkS-T execution
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModGeneMarkST::execute() {
    // Outfiles: file/path.faa, file/path.fnn
    // WARNING GeneMarkS-T assumes working directory as output right now
    std::string     temp_lst_file;      // Absolute path to output lst file from GeneMark (cwd)
    std::string     temp_faa_file;      // Absolute path to output FAA file from GeneMark (cwd)
    std::string     temp_fnn_file;      // Absolute path to output FNN file from GeneMark (cwd)
    std::string     temp_hmm_file;      // Absolute path to output HMM file from GeneMark (cwd)
    std::string     temp_gms_log_file;  // Absolute path to output log file from GeneMark (cwd)
    std::string     genemark_cmd;       // Terminal comand for GeneMarkS-T
    std::string     genemark_std_out;   // Standard output from GeneMark terminal run
    std::string     line;               // Temp line string
    int32           err_code;           // Terminal command error code
    TerminalData    terminalData;       // Terminal command structure (Defined in TerminalCommands.h)

    // Set temporary outputs (these will be moved to FINAL outpaths
    temp_faa_file     = PATHS(mpFileSystem->get_cur_dir(), mTranscriptomeFilename + FileSystem::EXT_FAA);
    temp_fnn_file     = PATHS(mpFileSystem->get_cur_dir(), mTranscriptomeFilename + FileSystem::EXT_FNN);
    temp_lst_file     = PATHS(mpFileSystem->get_cur_dir(), mTranscriptomeFilename + FileSystem::EXT_LST);
    temp_hmm_file     = PATHS(mpFileSystem->get_cur_dir(), GENEMARK_HMM_FILE);
    temp_gms_log_file = PATHS(mpFileSystem->get_cur_dir(), GENEMARK_LOG_FILE);

    genemark_cmd     = mExePath + " -faa -fnn " + mInHits;
    genemark_std_out = PATHS(mModOutDir, GENEMARK_STD_OUT);

    terminalData.command        = genemark_cmd;
    terminalData.print_files    = true;
    terminalData.base_std_path  = genemark_std_out;

    // !!!WARNING!!! GeneMarkS-T output is always in the CWD
    err_code = TC_execute_cmd(terminalData);
    if (err_code != 0 ) {
        throw ExceptionHandler("Error in running GeneMarkST at file located at: " +
                               mInHits + "\nGeneMarkST Error:\n" + terminalData.err_stream,
                               ERR_ENTAP_RUN_GENEMARK);
    }
    FS_dprint("Success!");

    // Ensure files successfully printed
    if (!mpFileSystem->file_exists(temp_faa_file)) {
        throw ExceptionHandler("Error unable to find GeneMarkS-T file located at: " + temp_faa_file,
                               ERR_ENTAP_RUN_GENEMARK);
    } else if (!mpFileSystem->file_exists(temp_fnn_file)) {
        throw ExceptionHandler("Error unable to find GeneMarkS-T file located at: " + temp_fnn_file,
                               ERR_ENTAP_RUN_GENEMARK);
    } else if (!mpFileSystem->file_exists(temp_lst_file)) {
        throw ExceptionHandler("Error unable to find GeneMarkS-T file located at: " + temp_lst_file,
                               ERR_ENTAP_RUN_GENEMARK);
    }

    FS_dprint("GeneMarkS-T files printed to:\nFAA File: " + temp_faa_file +
                                            "\nFNN File: " + temp_fnn_file+
                                            "\nLST File: " + temp_lst_file);

    // Format genemarks-t output (remove blank lines for user)
    // Format FNN file
    FS_dprint("Formatting GeneMarkST files...");
    try {
        std::ifstream in_file_fnn(temp_fnn_file);
        std::ofstream out_file_fnn(mFinalFnnPath);
        while(getline(in_file_fnn, line)) {
            if (!line.empty()) {
                out_file_fnn << line << '\n';
            }
        }
        in_file_fnn.close();
        out_file_fnn.close();

        // Format FAA file
        std::ifstream in_file_faa(temp_faa_file);
        std::ofstream out_file_faa(mFinalFaaPath);
        while(getline(in_file_faa, line)) {
            if (!line.empty()) {
                out_file_faa << line << '\n';
            }
        }
        in_file_faa.close();
        out_file_faa.close();
    } catch (std::exception &e) {
        throw ExceptionHandler("Error formatting GeneMarkS-T results\n" + std::string(e.what()),
                               ERR_ENTAP_RUN_GENEMARK_MOVE);
    }

    // Delete temporary files (ignore errors)
    mpFileSystem->delete_file(temp_fnn_file);
    mpFileSystem->delete_file(temp_faa_file);

    // Move other output files to module directory
    if (!mpFileSystem->rename_file(temp_lst_file, mFinalLstPath)) {
        throw ExceptionHandler("Error moving GeneMarkS-T results", ERR_ENTAP_RUN_GENEMARK_MOVE);
    }

    // move log file, ignore errors not needed for execution
    mpFileSystem->rename_file(temp_gms_log_file, mFinalGmstLogPath);

    // Does not always exist, not needed for final calculation so ignore errors
    mpFileSystem->rename_file(temp_hmm_file,mFinalHmmPath);
    FS_dprint("Success!");
}


/**
 * ======================================================================
 * Function void ModGeneMarkST::parse()
 *
 * Description          - Parses/calculates statistics of GeneMark output
 *                        and prints for user
 *                      - Calls Graphing Manager to graph results
 *                      - Prints genes (partial, complete, internal...) to different files
 *                      - Updates QueryData structure with new frame information (sequence, frame)
 *
 * Notes                - Fatal Error if any issue occurs
 *
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModGeneMarkST::parse() {
    // generate maps, query->sequence
    FS_dprint("Beginning to calculate Genemark statistics...");

    std::string                             out_removed_path;       // Absolute base path to removed sequences (no frame)
    std::string                             out_internal_path;      // Absolute base path to "internal" genes
    std::string                             out_complete_path;      // Absolute base path to "complete" genes
    std::string                             out_partial_path;       // Absolute base path to "partial" genes
    std::string                             figure_results_path;
    std::string                             figure_results_png;
    std::string                             figure_removed_path;
    std::string                             figure_removed_png;
    std::string                             min_removed_seq;        // Sequence ID of shortest removed sequence
    std::string                             min_kept_seq;           // Sequence ID of shortest kept sequence
    std::string                             max_removed_seq;        // Sequence ID of longest removed sequence
    std::string                             max_kept_seq;           // Sequence ID of longest kept sequence
    std::stringstream                       stat_output;            // Output string stream for statistics to Log File
    std::map<std::string, std::ofstream*>   file_map_faa;           // File map of FAA sequence and frame types
    std::map<std::string, std::ofstream*>   file_map_fnn;           // File map of FNN sequences ONLY (NO frame types)
    std::map<std::string, uint32>           count_map;              // Count of sequence types
    std::vector<uint16>                     all_kept_lengths;       // Total nucleotide length of kept sequences
    std::vector<uint16>                     all_lost_lengths;       // Total nucleotide length of removed sequences
    uint16                                  length;                 // Temp nucleotide length variable
    fp32                                    avg_selected;           // Average nucleotide length of kept sequences
    fp32                                    avg_lost;               // Average nucleotide length of removed sequences
    uint32                                  min_removed;            // bp count of smallest removed sequence
    uint32                                  min_selected;           // bp count of smallest kept sequence
    uint32                                  max_removed;            // bp count of largest removed sequence
    uint32                                  max_selected;           // bp count of largest kept sequence
    uint32                                  count_selected;         // Number of kept genes (frame was found)
    uint32                                  count_removed;          // Number of removed genes (frame was NOT found)
    uint64                                  total_removed_len;      // Total bp count of removed sequences
    uint64                                  total_kept_len;         // Total bp count of kept sequences
    std::pair<uint64, uint64>               kept_n;                 // N value of kept sequences
    GraphingData                            graphingStruct;         // Structure containing information for Graphing Manager

    // Ensure paths we need exist
    if (!mpFileSystem->file_exists(mFinalFaaPath)) {
        throw ExceptionHandler("Final GeneMarkST output not found at: " + mFinalFaaPath,
            ERR_ENTAP_RUN_GENEMARK_PARSE);
    } else if (!mpFileSystem->file_exists(mFinalLstPath)) {
        throw ExceptionHandler("Final GeneMarkST lst output not found at: " + mFinalLstPath,
                               ERR_ENTAP_RUN_GENEMARK_PARSE);
    } else if (!mpFileSystem->file_exists(mFinalFnnPath)) {
        throw ExceptionHandler("Final GeneMarkST lst output not found at: " + mFinalFnnPath,
                               ERR_ENTAP_RUN_GENEMARK_PARSE);
    }

    // Set up outpaths, directories are already created by super
    out_removed_path    = PATHS(mProcDir, FRAME_SELECTION_LOST);
    out_internal_path   = PATHS(mProcDir, FRAME_SELECTION_INTERNAL);
    out_complete_path   = PATHS(mProcDir, FRAME_SELECTION_COMPLTE);
    out_partial_path    = PATHS(mProcDir, FRAME_SELECTION_PARTIAL);
    figure_removed_path = PATHS(mFigureDir, GRAPH_TEXT_REF_COMPAR);
    figure_removed_png  = PATHS(mFigureDir, GRAPH_FILE_REF_COMPAR);
    figure_results_path = PATHS(mFigureDir, GRAPH_TEXT_FRAME_RESUTS);
    figure_results_png  = PATHS(mFigureDir, GRAPH_FILE_FRAME_RESUTS);

    // all nucleotide lengths
    min_removed=0xFFFFFFFF;
    min_selected=0xFFFFFFFF;
    max_removed=0;
    max_selected=0;
    total_removed_len=0;
    total_kept_len=0;

    // Sequence count
    count_selected=0;
    count_removed=0;

    try {
        // Parse protein file (.faa)
        std::map<std::string,FrameData> protein_map = genemark_parse_fasta(mFinalFaaPath);
        // Parse nucleotide file (.fnn) TODO cleanup/ move to query data just throwing in for now
        std::map<std::string,FrameData> nucleotide_map = genemark_parse_fasta(mFinalFnnPath);
        // Parse lst file to get info for each sequence (partial, internal...)
        genemark_parse_lst(mFinalLstPath,protein_map);

        std::ofstream file_figure_removed(figure_removed_path,std::ios::out | std::ios::app);
        std::ofstream file_figure_results(figure_results_path,std::ios::out | std::ios::app);

        file_figure_removed << "flag\tsequence length" << std::endl;    // First line placeholder, not used
        file_figure_results << "flag\tsequence length" << std::endl;


        // Initialize protein output streams (will be moved/changed)
        file_map_faa[FRAME_SELECTION_INTERNAL_FLAG] =
                new std::ofstream(out_internal_path+ FileSystem::EXT_FAA, std::ios::out | std::ios::app);
        file_map_faa[FRAME_SELECTION_COMPLETE_FLAG] =
                new std::ofstream(out_complete_path + FileSystem::EXT_FAA, std::ios::out | std::ios::app);
        file_map_faa[FRAME_SELECTION_FIVE_FLAG] =
                new std::ofstream(out_partial_path + FileSystem::EXT_FAA, std::ios::out | std::ios::app);
        file_map_faa[FRAME_SELECTION_THREE_FLAG] = file_map_faa[FRAME_SELECTION_FIVE_FLAG];

        // Initialize nucleotide output streams
        file_map_fnn[FRAME_SELECTION_LOST_FLAG] =
                new std::ofstream(out_removed_path + FileSystem::EXT_FNN, std::ios::out | std::ios::app);
        file_map_fnn[FRAME_SELECTION_INTERNAL_FLAG] =
                new std::ofstream(out_internal_path+ FileSystem::EXT_FNN, std::ios::out | std::ios::app);
        file_map_fnn[FRAME_SELECTION_COMPLETE_FLAG] =
                new std::ofstream(out_complete_path+ FileSystem::EXT_FNN, std::ios::out | std::ios::app);
        file_map_fnn[FRAME_SELECTION_FIVE_FLAG] =
                new std::ofstream(out_partial_path+ FileSystem::EXT_FNN, std::ios::out | std::ios::app);
        file_map_fnn[FRAME_SELECTION_THREE_FLAG] = file_map_fnn[FRAME_SELECTION_FIVE_FLAG];

        count_map ={
                {FRAME_SELECTION_INTERNAL_FLAG,0 },
                {FRAME_SELECTION_COMPLETE_FLAG,0 },
                {FRAME_SELECTION_FIVE_FLAG    ,0 },
                {FRAME_SELECTION_THREE_FLAG   ,0 },
        };

        for (auto& pair : *mpQueryData->get_sequences_ptr()) {
            std::map<std::string,FrameData>::iterator p_it = protein_map.find(pair.first);
            if (!pair.second->is_kept()) continue; // Skip seqs that were lost to expression analysis
            if (p_it != protein_map.end()) {
                // Kept sequence, either partial, complete, or internal
                count_selected++;
                pair.second->set_sequence_p(p_it->second.sequence); // Sets isprotein flag
                pair.second->setFrame(p_it->second.frame_type);

                // Update nucleotide sequence
                auto n_it = nucleotide_map.find(pair.first);
                if (n_it != nucleotide_map.end()) {
                    pair.second->set_sequence_n(n_it->second.sequence);
                }

                // Get Nucleotide sequence length for statistic purposes
                length = (uint16) pair.second->get_sequence_length();

                if (length < min_selected) {
                    min_selected = length;
                    min_kept_seq = pair.first;
                }
                if (length > max_selected) {
                    max_selected = length;
                    max_kept_seq = pair.first;
                }

                total_kept_len += length;            // Update total transcriptome length for statistics
                all_kept_lengths.push_back(length);  // Update individual lengths for statistics

                // Update figure
                file_figure_removed << GRAPH_KEPT_FLAG << '\t' << std::to_string(length) << std::endl;

                // Print nucleotide + protein sequences to files
                std::map<std::string, std::ofstream*>::iterator file_it = file_map_faa.find(p_it->second.frame_type);
                std::map<std::string, std::ofstream*>::iterator file_it_n = file_map_fnn.find(p_it->second.frame_type);
                if (file_it != file_map_faa.end() && file_it_n != file_map_faa.end()) {
                    *file_it->second << p_it->second.sequence << std::endl;
                    *file_it_n->second << n_it->second.sequence << std::endl;
                    count_map[p_it->second.frame_type]++;
                } else {
                    throw ExceptionHandler("Unknown frame flag found", ERR_ENTAP_RUN_GENEMARK_STATS);
                }

            } else {
                // Lost sequence from Frame Selection process
                count_removed++;
                pair.second->QUERY_FLAG_CLEAR(QuerySequence::QUERY_FRAME_KEPT);
                *file_map_fnn[FRAME_SELECTION_LOST_FLAG] << pair.second->get_sequence_n() << std::endl;

                length = (uint16) pair.second->get_sequence_length();  // Nucleotide sequence length for statistics

                if (length < min_removed) {
                    min_removed = length;
                    min_removed_seq = pair.first;
                }
                if (length > max_removed) {
                    max_removed_seq = pair.first;
                    max_removed = length;
                }
                file_figure_removed << GRAPH_REJECTED_FLAG << '\t' << std::to_string(length) << std::endl;
                all_lost_lengths.push_back(length);
                total_removed_len += length;
            }
        }

        // Cleanup/close protein files
        for(auto& pair : file_map_faa) {
            if (pair.first.compare(FRAME_SELECTION_THREE_FLAG)!=0){
                pair.second->close();
                delete pair.second;
                pair.second = 0;
            }
        }
        // Cleanup.close nucleotide files
        for(auto& pair : file_map_fnn) {
            if (pair.first.compare(FRAME_SELECTION_THREE_FLAG)!=0){
                pair.second->close();
                delete pair.second;
                pair.second = 0;
            }
        }

        // Ensure some sequences were kept and not all removed before we continue with printing stats, etc...
        if (count_selected <= MINIMUM_KEPT_SEQUENCES) {
            throw ExceptionHandler("Number of Kept sequences were under minimum after Frame Selection process, exiting...",
                ERR_ENTAP_RUN_GENEMARK_PARSE);
        }

        // Calculate and print stats
        FS_dprint("Beginning to calculate statistics...");
        avg_selected = (fp32)total_kept_len / count_selected;
        mpFileSystem->format_stat_stream(stat_output, "Frame Selected Transcripts (GeneMarkS-T)");
        stat_output <<
                    "Total sequences frame selected: "      << count_selected          <<
                    "\n\tTranslated protein sequences: "    << mFinalFaaPath              <<
                    "\nTotal sequences removed (no frame): "<< count_removed           <<
                    "\n\tFrame selected CDS removed: "      << out_removed_path        <<
                    "\nTotal of "                           <<
                    count_map[FRAME_SELECTION_FIVE_FLAG]    << " 5 prime partials and "<<
                    count_map[FRAME_SELECTION_THREE_FLAG]   << " 3 prime partials"     <<
                    "\n\tPartial CDS: "                     << out_partial_path        <<
                    "\nTotal of "                           <<
                    count_map[FRAME_SELECTION_COMPLETE_FLAG]<<" complete genes:\n\t" << out_complete_path<<
                    "\nTotal of "                           <<
                    count_map[FRAME_SELECTION_INTERNAL_FLAG]<<" internal genes:\n\t" << out_internal_path<<"\n\n";

        // Determine statistics for NEW reference transcriptome (following frame selection process)
        mpFileSystem->format_stat_stream(stat_output, "Frame Selection: New Reference Transcriptome Statistics");
        kept_n = mpQueryData->calculate_N_vals(all_kept_lengths,total_kept_len);
        stat_output <<
                    "\nTotal sequences: "      << count_selected <<
                    "\nTotal length of transcriptome(bp): "      << total_kept_len <<
                    "\nAverage length(bp): "   << avg_selected   <<
                    "\nn50: "                  << kept_n.first   <<
                    "\nn90: "                  << kept_n.second  <<
                    "\nLongest sequence(bp): " << max_selected   << " (" << max_kept_seq << ")" <<
                    "\nShortest sequence(bp): "<< min_selected   << " (" << min_kept_seq << ")";

        // Determine statistics for removed sequences following Frame Selection process
        if (count_removed > 0) {
            avg_lost     = (fp32)total_removed_len / count_removed;
            std::pair<uint64, uint64> removed_n =
                    mpQueryData->calculate_N_vals(all_lost_lengths,total_removed_len);
            stat_output <<
                        "\n\nRemoved Sequences (no frame):"       <<
                        "\nTotal sequences: "                     << count_removed    <<
                        "\nAverage sequence length(bp): "         << avg_lost         <<
                        "\nn50: "                                 << removed_n.first  <<
                        "\nn90: "                                 << removed_n.second <<
                        "\nLongest sequence(bp): "  << max_removed<< " (" << max_removed_seq << ")" <<
                        "\nShortest sequence(bp): " << min_removed<< " (" << min_removed_seq << ")" <<"\n";
        } else {
            stat_output << "WARNING: No sequences were removed from Frame Selection";
        }
        std::string stat_out_msg = stat_output.str();
        mpFileSystem->print_stats(stat_out_msg);
        FS_dprint("Success!");

        //---------------------- Figure handling ----------------------//
        FS_dprint("Beginning figure handling...");
        file_figure_results << GRAPH_REJECTED_FLAG           << '\t' << std::to_string(count_removed)   <<std::endl;
        file_figure_results << FRAME_SELECTION_FIVE_FLAG     << '\t' << std::to_string(count_map[FRAME_SELECTION_FIVE_FLAG]) <<std::endl;
        file_figure_results << FRAME_SELECTION_THREE_FLAG    << '\t' << std::to_string(count_map[FRAME_SELECTION_THREE_FLAG]) <<std::endl;
        file_figure_results << FRAME_SELECTION_COMPLETE_FLAG << '\t' << std::to_string(count_map[FRAME_SELECTION_COMPLETE_FLAG])   <<std::endl;
        file_figure_results << FRAME_SELECTION_INTERNAL_FLAG << '\t' << std::to_string(count_map[FRAME_SELECTION_INTERNAL_FLAG])   <<std::endl;

        graphingStruct.text_file_path = figure_results_path;
        graphingStruct.graph_title    = GRAPH_TITLE_FRAME_RESULTS;
        graphingStruct.fig_out_path   = figure_results_png;
        graphingStruct.software_flag  = GRAPH_FRAME_FLAG;
        graphingStruct.graph_type     = GRAPH_PIE_RESULTS_FLAG;
        mpGraphingManager->graph(graphingStruct);

        graphingStruct.text_file_path = figure_removed_path;
        graphingStruct.graph_title    = GRAPH_TITLE_REF_COMPAR;
        graphingStruct.fig_out_path   = figure_removed_png;
        graphingStruct.graph_type     = GRAPH_COMP_BOX_FLAG;
        mpGraphingManager->graph(graphingStruct);
        FS_dprint("Success!");
        //-----------------------------------------------------------//

    } catch (const ExceptionHandler &e) {throw e;}
    FS_dprint("Success! Parsing complete");
}


/**
 * ======================================================================
 * Function ModGeneMarkST::frame_map_t ModGeneMarkST::genemark_parse_protein
 *                                                  (std::string &protein)
 *
 * Description          - Analyzes fasta file to pull protein/nucleotide sequences that
 *                        will be added to QueryData
 * Notes                - None
 *
 * @param protein       - Absolute path to fasta (.faa or .fnn) file
 *
 * @return              - Structure containing protein sequence keyed to sequence ID
 *
 * =====================================================================
 */
ModGeneMarkST::frame_map_t ModGeneMarkST::genemark_parse_fasta(std::string &fasta) {
    FS_dprint("Parsing GeneMarkS-T FASTA file at: " + fasta);

    frame_map_t     protein_map;        // Map of new data from sequences  keyed to ID
    FrameData       protein_sequence;   // Data of translated sequences (sequence, frame type, etc...)
    std::string     line;               // Temp string
    std::string     sequence;           // FAA or FNN of translated sequence
    std::string     seq_id;             // ID of translated sequence
    uint16          first;
    uint16          second;

    // Filepath checked before
    std::ifstream in_file(fasta);
    while(true) {
        getline(in_file,line);
        if (line.empty() && !in_file.eof()) continue;
        if (line.find(FileSystem::FASTA_FLAG)==0 || in_file.eof()) {
            if (!seq_id.empty()) {
                if (in_file.eof()) {
                    sequence += line + "\n";
                }
                protein_sequence = {sequence, ""};
                protein_map.emplace(seq_id,protein_sequence);
            }
            if (in_file.eof()) break;
            // GeneMarkS-T adds a tab character within the sequence headers
            // Remove this to find the 'actual' sequence header we can reference against
            // our original transcriptome
            first    = (uint16) (line.find(FileSystem::FASTA_FLAG)+1);
            second   = (uint16) line.find('\t');
            if (first != std::string::npos || second != std::string::npos) {
                seq_id   = line.substr(first,line.find('\t')-first);
                sequence = FileSystem::FASTA_FLAG + seq_id + "\n";
            }
        } else {
            sequence += line + "\n";
        }
    }

    FS_dprint("Success! File parsed");
    in_file.close();
    return protein_map;
}


/**
 * ======================================================================
 * Function void ModGeneMarkST::genemark_parse_lst(std::string &lst_path,
 *                                                 frame_map_t &current_map)
 *
 * Description          - Parses .lst file produced from GeneMark run to get Frame types
 *                        (partial, complete, etc...)
 *                      - Will update the current_map with frame type information
 *
 * Notes                - None
 *
 * @param lst_path     - Path to .lst file produced from GeneMark
 * @param current_map  - Map being used to analyze GeneMark data
 *
 * @return              - None
 *
 * =====================================================================
 */
void ModGeneMarkST::genemark_parse_lst(std::string &lst_path, frame_map_t &current_map) {
    FS_dprint("Parsing GeneMarkS-T LST file at: " + lst_path);

    std::string     frame;
    std::string     line;
    std::string     seq_id;
    bool            prime_5;
    bool            prime_3;
    uint16          first;

    // Filepath checked before
    std::ifstream in_file(lst_path);
    while (getline(in_file,line)) {
        if (line.empty()) continue;
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        if (line.find("FASTA") == 0) {
            first = (uint16) (line.find(':') + 1);
            seq_id = line.substr(first);
        } else if (isdigit(line.at(0))) {
            prime_5 = line.find('<') != std::string::npos;
            prime_3 = line.find('>') != std::string::npos;
            if (prime_5 && prime_3) {
                frame = FRAME_SELECTION_INTERNAL_FLAG;
            } else if (!prime_5 && !prime_3) {
                frame = FRAME_SELECTION_COMPLETE_FLAG;
            } else if (prime_5) {
                frame = FRAME_SELECTION_FIVE_FLAG;
            } else frame = FRAME_SELECTION_THREE_FLAG;
            std::map<std::string, FrameData>::iterator it = current_map.find(seq_id);
            if (it != current_map.end()) {
                it->second.frame_type = frame;
            } else {
                throw ExceptionHandler("Sequence: " + seq_id + " not found in map during parsing of "
                                       "lst file at: " + mFinalLstPath,
                                       ERR_ENTAP_RUN_GENEMARK_PARSE);
            }
        }
    }
    FS_dprint("Success!");
}

/**
 * ======================================================================
 * Function void ModGeneMarkST::~ModGeneMarkST()
 *
 * Description          - Cleanup
 *
 * Notes                - Destructor
 *
 * @return              - None
 *
 * =====================================================================
 */
ModGeneMarkST::~ModGeneMarkST() {
    FS_dprint("Killing object - ModGeneMarkST");
}

/**
 * ======================================================================
 * Function void ModGeneMarkST::ModGeneMarkST()
 *
 * Description          - Initializes ModGeneMarkST object calling superclass AbstractFrame
 *                      - Sets final paths that GeneMarkS-T will create
 *
 * Notes                - Constructor
 *
 * @return              - None
 *
 * =====================================================================
 */
ModGeneMarkST::ModGeneMarkST(std::string &execution_stage_path, std::string &in_hits,
                             EntapDataPtrs &entap_data, std::string &exe) :
    AbstractFrame(execution_stage_path, in_hits, entap_data, "GeneMarkS-T", exe) {
    mTranscriptomeFilename = mpFileSystem->get_filename(in_hits, true);

    // Initialize FINAL output file paths
    // GenemarkS-T prints to CWD, they will be moved to these paths after execution
    mFinalFaaPath = PATHS(mModOutDir, mTranscriptomeFilename) + FileSystem::EXT_FAA;
    mFinalFnnPath = PATHS(mModOutDir, mTranscriptomeFilename) + FileSystem::EXT_FNN;
    mFinalLstPath = PATHS(mModOutDir, mTranscriptomeFilename + FileSystem::EXT_LST);
    mFinalGmstLogPath = PATHS(mModOutDir, GENEMARK_LOG_FILE);
    mFinalHmmPath = PATHS(mModOutDir, GENEMARK_HMM_FILE);
}

/**
 * ======================================================================
 * Function void ModGeneMarkST::get_final_faa()
 *
 * Description          - Returns absolute path to final faa (protein) file
 *                        created from GeneMarkS-T
 *
 * Notes                - None
 *
 * @return              - Absolute path to final faa (protein) file
 *
 * =====================================================================
 */
std::string ModGeneMarkST::get_final_faa() {
    return mFinalFaaPath;
}
