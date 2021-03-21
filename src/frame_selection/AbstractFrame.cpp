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

#include "AbstractFrame.h"
#include "../ExceptionHandler.h"

/**
 * ======================================================================
 * Function AbstractFrame(std::string &execution_stage_path, std::string &in_hits,
 *                        EntapDataPtrs &entap_data, std::string module_name,
                          std::string &exe)
 *
 * Description          - Constructor for Abstract frame selection class
 *                      - Initializes protected member variables for
 *                        expression modules
 *
 * Notes                - Constructor
 *
 * @param execution_stage_path - Absolute path to the output directory for this stage (Frame Selection)
 * @param in_hits              - Absolute path to input transcriptome
 * @param entap_data           - Pointers to necessary entap data for frame selection
 * @param module_name          - Name of this software module
 * @param exe                  - Execution method (i.e. executable)
 *
 * @return              - AbstractFrame object
 * ======================================================================
 */
AbstractFrame::AbstractFrame(std::string &execution_stage_path, std::string &in_hits,
                             EntapDataPtrs &entap_data, std::string module_name,
                             std::vector<ENTAP_HEADERS> &module_headers)
: EntapModule(execution_stage_path, in_hits, entap_data, module_name, module_headers) {

    mExecutionState = FRAME_SELECTION;
}

void AbstractFrame::frame_calculate_statistics() {
    FS_dprint("Beginning to calculate Frame Selection statistics...");

    std::string                             out_removed_path;       // Absolute base path to removed sequences (no frame)
    std::string                             out_internal_path;      // Absolute base path to "internal" genes
    std::string                             out_complete_path;      // Absolute base path to "complete" genes
    std::string                             out_partial_path;       // Absolute base path to "partial" genes
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

    // Set up outpaths, directories are already created by super
    out_removed_path    = PATHS(mProcDir, FRAME_SELECTION_FILENAME_LOST);
    out_internal_path   = PATHS(mProcDir, FRAME_SELECTION_FILENAME_INTERNAL);
    out_complete_path   = PATHS(mProcDir, FRAME_SELECTION_FILENAME_COMPLETE);
    out_partial_path    = PATHS(mProcDir, FRAME_SELECTION_FILENAME_PARTIAL);

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

        // Initialize graphing data
        GraphingManager::GraphingData graph_pie_results;
        graph_pie_results.x_axis_label = ENT_GRAPH_NULL;
        graph_pie_results.y_axis_label = "Sequence Length";
        graph_pie_results.text_file_path = PATHS(mFigureDir, GRAPH_TEXT_FRAME_RESUTS);
        graph_pie_results.fig_out_path   = PATHS(mFigureDir, GRAPH_FILE_FRAME_RESUTS);
        graph_pie_results.graph_title    = GRAPH_TITLE_FRAME_RESULTS;
        graph_pie_results.graph_type     = GraphingManager::ENT_GRAPH_PIE_CHART;

        GraphingManager::GraphingData graph_box_comparison;
        graph_box_comparison.x_axis_label = ENT_GRAPH_NULL;
        graph_box_comparison.y_axis_label = "Sequence Length";
        graph_box_comparison.text_file_path = PATHS(mFigureDir, GRAPH_TEXT_REF_COMPAR);
        graph_box_comparison.fig_out_path   = PATHS(mFigureDir, GRAPH_FILE_REF_COMPAR);
        graph_box_comparison.graph_title    = GRAPH_TITLE_REF_COMPAR;
        graph_box_comparison.graph_type     = GraphingManager::ENT_GRAPH_BOX_PLOT_VERTICAL;

        mpGraphingManager->initialize_graph_data(graph_pie_results);
        mpGraphingManager->initialize_graph_data(graph_box_comparison);


        for (auto& pair : *mpQueryData->get_sequences_ptr()) {
            if (!pair.second->is_kept()) continue; // Skip seqs that were lost to expression analysis

            if (pair.second->is_protein() && pair.second->is_kept_expression()) {
                // Kept sequence, either partial, complete, or internal
                count_selected++;

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
                mpGraphingManager->add_datapoint(graph_box_comparison.text_file_path, {GRAPH_KEPT_FLAG, std::to_string(length)});

                // Print nucleotide + protein sequences to files
                auto file_it = file_map_faa.find(pair.second->getFrame());
                auto file_it_n = file_map_fnn.find(pair.second->getFrame());
                if (file_it != file_map_faa.end() && file_it_n != file_map_faa.end()) {
                    *file_it->second << pair.second->get_sequence_p() << std::endl;
                    *file_it_n->second << pair.second->get_sequence_n() << std::endl;
                    count_map[pair.second->getFrame()]++;
                } else {
                    throw ExceptionHandler("Unknown frame flag found", ERR_ENTAP_RUN_GENEMARK_STATS);
                }

            } else if (!pair.second->is_protein() && pair.second->is_kept_expression()){
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
                mpGraphingManager->add_datapoint(graph_box_comparison.text_file_path, {GRAPH_REJECTED_FLAG, std::to_string(length)});
                all_lost_lengths.push_back(length);
                total_removed_len += length;
            } else {
                ;
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
        mpFileSystem->format_stat_stream(stat_output, "Frame Selected Transcripts (" + mModuleName + ")");
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
        mpGraphingManager->add_datapoint(graph_pie_results.text_file_path, {GRAPH_REJECTED_FLAG, std::to_string(count_removed)});
        mpGraphingManager->add_datapoint(graph_pie_results.text_file_path, {FRAME_SELECTION_FIVE_FLAG, std::to_string(count_map[FRAME_SELECTION_FIVE_FLAG])});
        mpGraphingManager->add_datapoint(graph_pie_results.text_file_path, {FRAME_SELECTION_THREE_FLAG, std::to_string(count_map[FRAME_SELECTION_THREE_FLAG])});
        mpGraphingManager->add_datapoint(graph_pie_results.text_file_path, {FRAME_SELECTION_COMPLETE_FLAG, std::to_string(count_map[FRAME_SELECTION_COMPLETE_FLAG])});
        mpGraphingManager->add_datapoint(graph_pie_results.text_file_path, {FRAME_SELECTION_INTERNAL_FLAG, std::to_string(count_map[FRAME_SELECTION_INTERNAL_FLAG])});

        mpGraphingManager->graph_data(graph_pie_results.text_file_path);
        mpGraphingManager->graph_data(graph_box_comparison.text_file_path);
        //------------------------------------------------------------//

        FS_dprint("Success! Frame Selection statistics completed");

    } catch (const std::exception &e) {
        throw ExceptionHandler("ERROR Unable to calculate Frame Selection statistics: " + std::string(e.what()),
                               ERR_ENTAP_RUN_GENEMARK_PARSE);
    }
}


/**
 * ======================================================================
 * Function void AbstractFrames::get_final_faa()
 *
 * Description          - Returns absolute path to final faa (protein) file
 *                        created from Frame Selection software
 *
 * Notes                - None
 *
 * @return              - Absolute path to final faa (protein) file
 *
 * =====================================================================
 */
std::string AbstractFrame::get_final_faa() {
    return mFinalFaaPath;
}

void AbstractFrame::set_success_flags() {
    mpQueryData->set_is_protein_data(true);
    mpQueryData->set_is_success_frame_selection(true);
}

