//
// Created by harta on 8/11/17.
//

#ifndef ENTAP_MODGENEMARKST_H
#define ENTAP_MODGENEMARKST_H


#include "AbstractFrame.h"

class ModGeneMarkST : public AbstractFrame{

struct frame_seq {
    unsigned long length;
    std::string sequence;
    std::string frame_type;
};
typedef std::map<std::string,ModGeneMarkST::frame_seq> frame_map_t;

public:
    ModGeneMarkST(std::string &exe, std::string &out, std::string &in, std::string &proc,
                  std::string &fig, std::string &frame, GraphingManager *graphing) :
            AbstractFrame(exe, out, in, proc, fig, frame, graphing){}

    virtual std::pair<bool, std::string> verify_files() override ;
    virtual std::string execute(std::map<std::string, QuerySequence>&) override ;
    virtual void parse(std::map<std::string, QuerySequence>&) override ;

private:

    const std::string GRAPH_TITLE_FRAME_RESULTS     = "Frame_Selection_ORFs";
    const std::string GRAPH_FILE_FRAME_RESUTS       = "frame_results_pie.png";
    const std::string GRAPH_TEXT_FRAME_RESUTS       = "frame_results_pie.txt";
    const std::string GRAPH_TITLE_REF_COMPAR        = "Frame_Selected_Sequences";
    const std::string GRAPH_FILE_REF_COMPAR         = "rejected_comparison_box.png";
    const std::string GRAPH_TEXT_REF_COMPAR         = "rejected_comparison_box.txt";
    const std::string GRAPH_REJECTED_FLAG           = "Removed";
    const std::string GRAPH_KEPT_FLAG               = "Selected";

    const unsigned char GRAPH_FRAME_FLAG            = 1;    // Ensure these match with entap_graphing.py
    const unsigned char GRAPH_PIE_RESULTS_FLAG      = 1;
    const unsigned char GRAPH_COMP_BOX_FLAG         = 2;

    const std::string GENEMARK_LOG_FILE             = "gms.log";
    const std::string GENEMARK_HMM_FILE             = "GeneMark_hmm.mod";
    const std::string GENEMARK_STD_OUT              = "genemark_run";
    const std::string FRAME_SELECTION_PARTIAL       = "partial_genes.fasta";
    const std::string FRAME_SELECTION_COMPLTE       = "complete_genes.fasta";
    const std::string FRAME_SELECTION_INTERNAL      = "internal_genes.fasta";
    const std::string FRAME_SELECTION_LOST          = "sequences_removed.fasta";
    const std::string FRAME_SELECTION_LOST_FLAG     = "lost";
    const std::string FRAME_SELECTION_FIVE_FLAG     = "Partial 5 Prime";
    const std::string FRAME_SELECTION_THREE_FLAG    = "Partial 3 Prime";
    const std::string FRAME_SELECTION_COMPLETE_FLAG = "Complete";
    const std::string FRAME_SELECTION_INTERNAL_FLAG = "Internal";

    std::string _final_out;
    std::string _final_lst;


    frame_map_t genemark_parse_protein(std::string&);
    void genemark_parse_lst(std::string &, std::map<std::string, frame_seq>&);
};


#endif //ENTAP_MODGENEMARKST_H
