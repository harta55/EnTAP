//
// Created by harta on 5/22/17.
//

#include <boost/filesystem.hpp>
#include <csv.h>
#include "Ontology.h"
#include "EntapConsts.h"
#include "EntapInit.h"
#include "ExceptionHandler.h"
#include "EntapExecute.h"

namespace boostFS = boost::filesystem;

Ontology::Ontology(int thread, bool overwrite, std::string egg_exe, std::string outpath,
                   std::string entap_exe) {
    _eggnog_exe = egg_exe;
    _threads = thread;
    _is_overwrite = overwrite;
    _entap_exe = entap_exe;
    _outpath = outpath;
}


void Ontology::execute(short software, query_map_struct &SEQUENCES, std::string input,
                       std::string no_hit) {
    _new_input = input;
    _input_no_hits = no_hit;
    try {
        switch(software) {
            case 0:
                run_eggnog(SEQUENCES);
                break;
            default:
                run_eggnog(SEQUENCES);
                break;
        }
    } catch (ExceptionHandler &e) {throw e;}

}

/**
 *
 * @param SEQUENCES
 * @param out - <egg results of sequences that hit databases, results of no hits>
 */
void Ontology::parse_results_eggnog(query_map_struct& SEQUENCES, std::pair<std::string,std::string>& out) {
//    for (std::string data : out) {
//        entapInit::print_msg("Eggnog file located at " + data + " being filtered");
//        std::string msg = ENTAP_STATS::SOFTWARE_BREAK + "Ontology - Eggnog\n" +
//                          ENTAP_STATS::SOFTWARE_BREAK;
//        entapExecute::print_statistics(msg, _outpath);
//
//        io::CSVReader<ENTAP_EXECUTE::EGGNOG_COL_NUM, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(data);
//        // todo have columns from input file, in_read_header for versatility
//        std::string qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
//                best_og, cog_cat, eggnog_annot;
//        unsigned long count_go_hits, count_kegg_hits, count_total_go, count_total_kegg, count_no_hits;
//
//        std::map<std::string, QuerySequence> database_map;
//
//        while (in.read_row(qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
//                           best_og, cog_cat, eggnog_annot)) {
//        }
//
//    }
}


void Ontology::run_eggnog(query_map_struct &SEQUENCES) {

    std::string eggnog_out_dir = _outpath + ENTAP_EXECUTE::ONTOLOGY_OUT_PATH;
    boostFS::create_directories(eggnog_out_dir);
    std::string annotation_flag = eggnog_out_dir + "annotation_results";
    std::string annotation_std = eggnog_out_dir + "annotation_std";
    // TODO eggnog overwrite


    std::unordered_map<std::string,std::string> eggnog_command_map = {
            {"-i",_new_input},
            {"--output",annotation_flag},
            {"--cpu",std::to_string(_threads)},
            {"-m", "diamond"}
    };

    std::string eggnog_command = "python " + _eggnog_exe;
    for (auto &pair : eggnog_command_map) {
        eggnog_command += pair.first + " " + pair.second;
    }
    entapInit::print_msg("\nExecuting eggnog mapper against protein sequences that hit databases...\n"
                         + eggnog_command);
    if (entapInit::execute_cmd(eggnog_command, annotation_std) !=0) {
        throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
    }
    entapInit::print_msg("Success! Results written to: ");
    eggnog_command = "python " + _eggnog_exe;

    std::pair<std::string,std::string> out;
    parse_results_eggnog(SEQUENCES, out);
}

