//
// Created by harta on 5/7/17.
//

#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <boost/regex.hpp>
#include <iomanip>
#include "FrameSelection.h"
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapConfig.h"
#include "EntapGlobals.h"
#include "frame_selection/ModGeneMarkST.h"

namespace boostFS = boost::filesystem;

// can accept version/other for dependency injection
FrameSelection::FrameSelection(std::string &input, std::string &out,
                               boost::program_options::variables_map &user_flags,
                               GraphingManager *graphingManager) {
    print_debug("Spawn object - FrameSelection");
    _graphingManager = graphingManager;
    _exe_path        = GENEMARK_EXE;
    _outpath         = out;
    _inpath          = input;
    _overwrite       = (bool) user_flags.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _software_flag   = 0;
    _processed_path  = (boostFS::path(out) / boostFS::path(FRAME_SELECTION_OUT_DIR) /
                        boostFS::path(PROCESSED_DIR)).string();
    _figure_path     = (boostFS::path(out) / boostFS::path(FRAME_SELECTION_OUT_DIR) /
                        boostFS::path(FIGURE_DIR)).string();
    _frame_outpath   = (boostFS::path(out) / boostFS::path(FRAME_SELECTION_OUT_DIR)).string();
    SOFTWARE = static_cast<FrameSoftware >(_software_flag);
}


std::string FrameSelection::execute(std::string input, std::map<std::string,QuerySequence> &SEQUENCES) {

    std::string output;
    std::pair<bool, std::string> verify_pair;

    _inpath = input;
    if (_overwrite) boostFS::remove_all(_frame_outpath);
    boostFS::create_directories(_frame_outpath);
    boostFS::create_directories(_processed_path);
    boostFS::create_directories(_figure_path);
    try {
        std::unique_ptr<AbstractFrame> ptr = spawn_object();
        verify_pair = ptr->verify_files();
        if (!verify_pair.first) {
            output = ptr->execute(SEQUENCES);
        } else output = verify_pair.second;
        ptr->parse(SEQUENCES);
        ptr.release();
        return output;
    } catch (const ExceptionHandler &e) {throw e;}
}


std::unique_ptr<AbstractFrame> FrameSelection::spawn_object() {
    // Handle any special conditions for each software

    switch (SOFTWARE) {
        case GENEMARKST:
            return std::unique_ptr<AbstractFrame>(new ModGeneMarkST(
                    _exe_path, _outpath, _inpath, _processed_path, _figure_path,
                    _frame_outpath, _graphingManager
            ));
        default:
            return std::unique_ptr<AbstractFrame>(new ModGeneMarkST(
                    _exe_path, _outpath, _inpath, _processed_path, _figure_path,
                    _frame_outpath, _graphingManager
            ));
    }
}

