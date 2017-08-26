//
// Created by harta on 5/22/17.
//

#include <boost/filesystem.hpp>
#include <csv.h>
#include <boost/archive/binary_iarchive.hpp>
#include "Ontology.h"
#include "EntapConfig.h"
#include <boost/serialization/map.hpp>
#include <iomanip>
#include "ExceptionHandler.h"
#include "EntapExecute.h"
#include "SimilaritySearch.h"
#include "DatabaseHelper.h"
#include "EntapGlobals.h"
#include "ontology/AbstractOntology.h"
#include "ontology/ModEggnog.h"

Ontology::Ontology(int thread, std::string outpath, std::string input,
                   boost::program_options::variables_map &user_input, GraphingManager* graphing) {
    print_debug("Spawn object - Ontology");
    _ontology_exe = EGG_EMAPPER_EXE;
    _threads = thread;
    _outpath = outpath;
    _new_input = input;
    _is_overwrite = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _software_flag = user_input[ENTAP_CONFIG::INPUT_FLAG_ONTOLOGY].as<short>();
    _go_levels = user_input[ENTAP_CONFIG::INPUT_FLAG_GO_LEVELS].as<std::vector<short>>();
    std::vector<std::string> _interpro_databases =
            user_input[ENTAP_CONFIG::INPUT_FLAG_INTERPRO].as<std::vector<std::string>>();
    _ontology_dir = (boostFS::path(outpath) / boostFS::path(ONTOLOGY_OUT_PATH)).string();
    _processed_dir = (boostFS::path(_ontology_dir) / boostFS::path(PROCESSED_OUT_DIR)).string();
    _figure_dir = (boostFS::path(_processed_dir) / boostFS::path(FIGURE_DIR)).string();
    _eggnog_db_path = EGG_SQL_DB_PATH;
    _graphingManager = graphing;
    SOFTWARE = static_cast<OntologySoftware>(_software_flag);
}


void Ontology::execute(query_map_struct &SEQUENCES, std::string input,std::string no_hit) {

    std::pair<bool,std::string> verify_pair;

    _new_input = input;
    _input_no_hits = no_hit;

    if (_is_overwrite) boostFS::remove_all(_ontology_dir);
    boostFS::create_directories(_ontology_dir);
    init_headers();
    try {
        std::unique_ptr<AbstractOntology> ptr = spawn_object();
        ptr->set_data(_go_levels,_eggnog_db_path,_threads);
        verify_pair = ptr->verify_files();
        if (!verify_pair.first) ptr->execute(SEQUENCES);
        ptr->parse(SEQUENCES);
        print_eggnog(SEQUENCES);

        // TODO move printing to manager
    } catch (ExceptionHandler &e) {throw e;}
}


std::unique_ptr<AbstractOntology> Ontology::spawn_object() {
    // Handle any special conditions for each software

    // Each will have separate outpaths implement soon...

    std::string proc_path;
    std::string exe_path;
    std::string fig_path;
    std::string out_path;

    switch (SOFTWARE) {
        case EGGNOG:
            return std::unique_ptr<AbstractOntology>(new ModEggnog(
                    _ontology_exe, _outpath, _new_input, _input_no_hits,
                    _processed_dir, _figure_dir, _ontology_dir, _graphingManager
            ));
        case INTERPRO:
            break;
        default:
            return std::unique_ptr<AbstractOntology>(new ModEggnog(
                    _ontology_exe, _outpath, _new_input, _input_no_hits,
                    _processed_dir, _figure_dir, _ontology_dir, _graphingManager
            ));
    }
}



void Ontology::print_eggnog(query_map_struct &SEQUENCES) {
    print_debug("Beginning to print final results...");
    std::map<short, std::ofstream*> file_map;
    std::string file_name;
    std::string outpath;
    for (short lvl : _go_levels) {
        file_name = "final_annotations_lvl" + std::to_string(lvl) + ".tsv";
        outpath = (boostFS::path(_outpath) / boostFS::path(file_name)).string();
        boostFS::remove(outpath);
        file_map[lvl] =
                new std::ofstream(outpath, std::ios::out | std::ios::app);
        for (const std::string *header : _HEADERS) {
            *file_map[lvl] << *header << '\t';
        }
        *file_map[lvl] << std::endl;
    }
    for (auto &pair : SEQUENCES) {
        for (short i : _go_levels) {
            *file_map[i]<< pair.second.print_tsv(_software_flag,_HEADERS,i)<<std::endl;
        }
    }
    for(auto& pair : file_map) {
        pair.second->close();
        delete pair.second;
    }
    print_debug("Success!");
}


void Ontology::init_headers() {
    switch (_software_flag) {
        case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
            _HEADERS = {
                    &ENTAP_EXECUTE::HEADER_QUERY,
                    &ENTAP_EXECUTE::HEADER_SUBJECT,
                    &ENTAP_EXECUTE::HEADER_PERCENT,
                    &ENTAP_EXECUTE::HEADER_ALIGN_LEN,
                    &ENTAP_EXECUTE::HEADER_MISMATCH,
                    &ENTAP_EXECUTE::HEADER_GAP_OPEN,
                    &ENTAP_EXECUTE::HEADER_QUERY_S,
                    &ENTAP_EXECUTE::HEADER_QUERY_E,
                    &ENTAP_EXECUTE::HEADER_SUBJ_S,
                    &ENTAP_EXECUTE::HEADER_SUBJ_E,
                    &ENTAP_EXECUTE::HEADER_E_VAL,
                    &ENTAP_EXECUTE::HEADER_COVERAGE,
                    &ENTAP_EXECUTE::HEADER_TITLE,
                    &ENTAP_EXECUTE::HEADER_SPECIES,
                    &ENTAP_EXECUTE::HEADER_DATABASE,
                    &ENTAP_EXECUTE::HEADER_FRAME,
                    &ENTAP_EXECUTE::HEADER_CONTAM,
                    &ENTAP_EXECUTE::HEADER_INFORM,
                    &ENTAP_EXECUTE::HEADER_SEED_ORTH,
                    &ENTAP_EXECUTE::HEADER_SEED_EVAL,
                    &ENTAP_EXECUTE::HEADER_SEED_SCORE,
                    &ENTAP_EXECUTE::HEADER_PRED_GENE,
                    &ENTAP_EXECUTE::HEADER_TAX_SCOPE,
                    &ENTAP_EXECUTE::HEADER_EGG_OGS,
                    &ENTAP_EXECUTE::HEADER_EGG_DESC,
                    &ENTAP_EXECUTE::HEADER_EGG_KEGG,
                    &ENTAP_EXECUTE::HEADER_EGG_PROTEIN,
                    &ENTAP_EXECUTE::HEADER_EGG_GO_BIO,
                    &ENTAP_EXECUTE::HEADER_EGG_GO_CELL,
                    &ENTAP_EXECUTE::HEADER_EGG_GO_MOLE
            };
            break;
        case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
//            _HEADERS = {
//                "Protein Database"      ,
//                "Protein ID"            ,
//                "Protein Term"          ,
//                "E_value"               ,
//                "InterPro ID"           ,
//                "InterPro Term"         ,
//                "GO Biological"         ,
//                "GO Cellular"           ,
//                "GO Molecular"          ,
//                "Pathway Terms"
//            };
            break;
        default:
            break;
    }
}

// TODO combine with eggnog
void Ontology::print_interpro(query_map_struct &SEQUENCES) {

//    std::string final_annotations;
//
//    final_annotations = (boostFS::path(_outpath) / boostFS::path("final_annotations.tsv")).string();
//    print_header(final_annotations);
//    std::ofstream file(final_annotations, std::ios::out | std::ios::app);
//    for (auto &pair : SEQUENCES) {
//        file<< pair.second.print_final_results(_software_flag,_HEADERS,0)<<std::endl;
//    }
//    file.close();
}

void Ontology::print_header(std::string file) {
    std::ofstream ofstream(file, std::ios::out | std::ios::app);
    for (const std::string *val : _HEADERS) ofstream << *val << '\t';
    ofstream<<std::endl;
    ofstream.close();
}
