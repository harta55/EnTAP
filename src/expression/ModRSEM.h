//
// Created by harta on 8/11/17.
//

#ifndef ENTAP_MODRSEM_H
#define ENTAP_MODRSEM_H


#include "AbstractExpression.h"

class ModRSEM : public AbstractExpression{

public:
    ModRSEM(std::string &exe, std::string &out, std::string &in, std::string &proc,
            std::string &fig, std::string &exp, std::string &align,GraphingManager *graphing) :
    AbstractExpression(exe, out, in, proc, fig, exp, align,graphing){}

    virtual std::pair<bool, std::string> verify_files() override ;
    virtual void execute(std::map<std::string, QuerySequence>&) override ;
    virtual std::string filter(std::map<std::string, QuerySequence>&) override ;
    virtual void set_data(int, float, bool) override    ;

private:

    const std::string GRAPH_TXT_BOX_PLOT    = "comparison_box.txt";
    const std::string GRAPH_PNG_BOX_PLOT    = "comparison_box.png";
    const std::string GRAPH_TITLE_BOX_PLOT  = "Expression_Analysis";
    const std::string GRAPH_REJECTED_FLAG   = "Removed";
    const std::string GRAPH_KEPT_FLAG       = "Selected";
    const float REJECTED_ERROR_CUTOFF       = 75.0;

    const unsigned char GRAPH_EXPRESSION_FLAG = 2;
    const unsigned char GRAPH_BOX_FLAG        = 1;
    static constexpr int RSEM_COL_NUM = 7;

    std::string _filename;
    std::string _rsem_out;
    std::string _exp_out;
    int         _threads;
    float       _fpkm;
    bool        _ispaired;

    bool rsem_validate_file(std::string);
    bool rsem_conv_to_bam(std::string);
    bool is_file_empty(std::string);
};


#endif //ENTAP_MODRSEM_H
