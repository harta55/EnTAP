/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/

//
//#ifndef ENTAP_MODINTERPRO_H
//#define ENTAP_MODINTERPRO_H
//
//
//#include "AbstractOntology.h"
//
//class ModInterpro : public AbstractOntology{
//
//    struct interpro_struct {
//        double _eval;
//        std::map<std::string,std::string> _results;
//        std::map<std::string,std::vector<std::string>> _go_map;
//    };
//
//public:
//    ModInterpro(std::string &exe, std::string &entap,std::string &out, std::string &in,
//            std::string &in_no_hits,std::string &proc,
//            std::string &fig, std::string &ont, GraphingManager *graphing) :
//    AbstractOntology(exe, entap, out, in, in_no_hits,proc, fig, ont, graphing){}
//
//    virtual std::pair<bool, std::string> verify_files() override ;
//    virtual void execute(std::map<std::string, QuerySequence>&) override ;
//    virtual void parse(std::map<std::string, QuerySequence>&) override ;
//    virtual void set_data(std::vector<short>&, std::string&, int) override ;
//
//
//private:
//    static constexpr short INTERPRO_COL_NUM = 15;
//
//};
//
//
//#endif //ENTAP_MODINTERPRO_H
