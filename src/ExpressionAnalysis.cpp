//
// Created by harta on 5/7/17.
//

#include <boost/filesystem.hpp>
#include <csv.h>
#include <boost/regex.hpp>
#include <iomanip>
#include "ExpressionAnalysis.h"
#include "ExceptionHandler.h"
#include "EntapConfig.h"
#include "EntapGlobals.h"
#include "EntapExecute.h"
#include "expression/ModRSEM.h"

namespace boostFS = boost::filesystem;

ExpressionAnalysis::ExpressionAnalysis(std::string &input,int t, std::string &exe, std::string &out
    , boost::program_options::variables_map& user_flags, GraphingManager *graph) {
    print_debug("Spawn object - ExpressionAnalysis");
    _inpath = input;
    _threads = t;
    _exepath = exe;
    _outpath = out;
    _software_flag = 0;
    _overwrite = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _ispaired = (bool) user_flags.count("paired-end");
    if (user_flags.count(ENTAP_CONFIG::INPUT_FLAG_ALIGN)) {
        _alignpath = user_flags[ENTAP_CONFIG::INPUT_FLAG_ALIGN].as<std::string>();
    }
    _fpkm = user_flags[ENTAP_CONFIG::INPUT_FLAG_FPKM].as<float>();
    _rsem_dir = (boostFS::path(out) / boostFS::path(RSEM_OUT_DIR)).string();
    _proc_dir = (boostFS::path(_rsem_dir) / boostFS::path(RSEM_PROCESSED_DIR)).string();
    _figure_dir = (boostFS::path(_proc_dir) / boostFS::path(RSEM_FIGURE_DIR)).string();
    _graphingManager = graph;
    SOFTWARE = static_cast<ExpressionSoftware>(_software_flag);
}

std::string ExpressionAnalysis::execute(std::string input,
                                        std::map<std::string, QuerySequence>& MAP) {
    std::string output;
    std::pair<bool, std::string> verify_pair;
    std::unique_ptr<AbstractExpression> ptr;

    _inpath = input;
    if (_overwrite) boostFS::remove_all(_rsem_dir);
    boostFS::create_directories(_rsem_dir);
    boostFS::create_directories(_figure_dir);
    boostFS::create_directories(_proc_dir);
    try {
        ptr = spawn_object();
        ptr->set_data(_threads, _fpkm, _ispaired);
        verify_pair = ptr->verify_files();
        if (!verify_pair.first) ptr->execute(MAP);
        output = ptr->filter(MAP);
        ptr.release();
    } catch (const ExceptionHandler &e) {throw e;}
    return output;
}

std::unique_ptr<AbstractExpression> ExpressionAnalysis::spawn_object() {
    switch (SOFTWARE) {
        case RSEM:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _exepath, _outpath, _inpath, _proc_dir, _figure_dir,
                    _rsem_dir, _alignpath, _graphingManager
            ));
        default:
            return std::unique_ptr<AbstractExpression>(new ModRSEM(
                    _exepath, _outpath, _inpath, _proc_dir, _figure_dir,
                    _rsem_dir, _alignpath, _graphingManager
            ));
    }
}


