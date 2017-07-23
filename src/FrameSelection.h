//
// Created by harta on 5/7/17.
//

#ifndef ENTAP_FRAMESELECTION_H
#define ENTAP_FRAMESELECTION_H
#include <iostream>
#include <map>
#include <boost/program_options/variables_map.hpp>
#include "QuerySequence.h"
#include "GraphingManager.h"


class FrameSelection {
    struct frame_seq {
        unsigned long length;
        std::string sequence;
        std::string frame_type;
    };
    typedef std::map<std::string,FrameSelection::frame_seq> frame_map_t;

    public:
        std::string execute(std::string,std::map<std::string,QuerySequence>&);
        FrameSelection(std::string&, std::string&, std::string&,
                       boost::program_options::variables_map &, GraphingManager*);

    private:

        const std::string GRAPH_TITLE_FRAME_RESULTS     = "Frame_Selection_Results";
        const std::string GRAPH_FILE_FRAME_RESUTS       = "frame_results_pie.png";
        const std::string GRAPH_TEXT_FRAME_RESUTS       = "frame_results_pie.txt";
        const std::string GRAPH_TITLE_REF_COMPAR        = "Rejected_Sequence_Comparison";
        const std::string GRAPH_FILE_REF_COMPAR         = "rejected_comparison_box.png";
        const std::string GRAPH_TEXT_REF_COMPAR         = "rejected_comparison_box.txt";
        const std::string GRAPH_REJECTED_FLAG           = "Rejected";
        const std::string GRAPH_KEPT_FLAG               = "Kept";
        const std::string FRAME_SELECTION_OUT_DIR       = "frame_selection/";
        const std::string FRAME_SELECTION_PARTIAL       = "partial_genes.fasta";
        const std::string FRAME_SELECTION_COMPLTE       = "complete_genes.fasta";
        const std::string FRAME_SELECTION_INTERNAL      = "internal_genes.fasta";
        const std::string FRAME_SELECTION_LOST          = "sequences_lost.fasta";
        const std::string FRAME_SELECTION_LOST_FLAG     = "lost";
        const std::string FRAME_SELECTION_FIVE_FLAG     = "Partial 5 Prime";
        const std::string FRAME_SELECTION_THREE_FLAG    = "Partial 3 Prime";
        const std::string FRAME_SELECTION_COMPLETE_FLAG = "Complete";
        const std::string FRAME_SELECTION_INTERNAL_FLAG = "Internal";
        const std::string GENEMARK_LOG_FILE             = "gms.log";
        const std::string GENEMARK_HMM_FILE             = "GeneMark_hmm.mod";
        const std::string GENEMARK_STD_OUT              = "genemark_run";

        const unsigned char GRAPH_FRAME_FLAG            = 1;    // Ensure these match with entap_graphing.py
        const unsigned char GRAPH_PIE_RESULTS_FLAG      = 1;
        const unsigned char GRAPH_COMP_BOX_FLAG         = 2;

        std::string      _frame_outpath;
        std::string      _processed_path;
        std::string      _figure_path;
        std::string      _exe_path;
        std::string      _inpath;
        std::string      _outpath;
        bool             _overwrite;
        short            _software_flag;
        GraphingManager  *_graphingManager;

        std::string genemarkst(std::map<std::string,QuerySequence> &);
        void genemarkStats(std::string&,std::string&,std::map<std::string,QuerySequence> &);
        frame_map_t genemark_parse_protein(std::string&);
        void genemark_parse_lst(std::string &, std::map<std::string, frame_seq>&);
};


#endif //ENTAP_FRAMESELECTION_H
