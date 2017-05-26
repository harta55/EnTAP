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
    entapInit::print_msg("Eggnog file located at " + out.first + " being filtered");
    std::string msg = ENTAP_STATS::SOFTWARE_BREAK + "Ontology - Eggnog\n" +
                      ENTAP_STATS::SOFTWARE_BREAK;
    entapExecute::print_statistics(msg, _outpath);
    std::string qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
            best_og, cog_cat, eggnog_annot;

    io::CSVReader<ENTAP_EXECUTE::EGGNOG_COL_NUM, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(out.first);
    entapInit::print_msg("Parsing results from sequences that hit databases");

    std::map<std::string, QuerySequence> database_map;
    in.next_line();in.next_line();in.next_line();
    while (in.read_row(qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                       best_og, cog_cat, eggnog_annot)) {
        SEQUENCES.at(qseqid).set_eggnog_results(seed_ortho,seed_e,seed_score,predicted_gene,go_terms,
            kegg,tax_scope,ogs);
    }
    // TODO stats
    print_eggnog(SEQUENCES);
}


void Ontology::run_eggnog(query_map_struct &SEQUENCES) {

    std::string eggnog_out_dir = _outpath + ENTAP_EXECUTE::ONTOLOGY_OUT_PATH;
    boostFS::create_directories(eggnog_out_dir);
//    std::string annotation_base_flag = eggnog_out_dir + "annotation_results";
//    std::string annotation_no_flag = eggnog_out_dir + "annotation_results_no_hits";
        std::string annotation_base_flag = "annotation_results";
    std::string annotation_no_flag = "annotation_results_no_hits";
    std::string annotation_std = eggnog_out_dir + "annotation_std";
    std::string eggnog_command = "python " + _eggnog_exe;
    std::pair<std::string,std::string> out;

    // TODO eggnog overwrite
    if (!entapInit::file_exists(_new_input) || !entapInit::file_exists(_input_no_hits)) {
        throw ExceptionHandler("File located at: " + _new_input +
                                       " or: " + _input_no_hits + "not found",ENTAP_ERR::E_RUN_EGGNOG);
    }
    std::unordered_map<std::string,std::string> eggnog_command_map = {
            {"-i",_new_input},
            {"--output",annotation_base_flag},
            {"--cpu",std::to_string(_threads)},
            {"-m", "diamond"}
    };

    for (auto &pair : eggnog_command_map) {
        eggnog_command += pair.first + " " + pair.second;
    }
    entapInit::print_msg("\nExecuting eggnog mapper against protein sequences that hit databases...\n"
                         + eggnog_command);
    if (entapInit::execute_cmd(eggnog_command, annotation_std) !=0) {
        throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
    }
    entapInit::print_msg("Success! Results written to: " + annotation_base_flag);
    eggnog_command_map["-i"] = _input_no_hits;
    eggnog_command_map["--output"] = annotation_no_flag;
    entapInit::print_msg("\nExecuting eggnog mapper against protein sequences that did not hit databases...\n"
                         + eggnog_command);
    if (entapInit::execute_cmd(eggnog_command, annotation_std) !=0) {
        throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
    }
    out.first = annotation_base_flag + ".emapper.annotations";
    out.second = annotation_no_flag + ".emapper.annotations";

    parse_results_eggnog(SEQUENCES, out);
}

void Ontology::print_eggnog(Ontology::query_map_struct &SEQUENCES) {
    std::string final_annotations = _outpath + "final_annotations.tsv";
    std::ofstream file(final_annotations, std::ios::out | std::ios::app);
    file <<
         "Query Seq\t"
                 "Subject Seq\t"
                 "Percent Identical\t"
                 "Alignment Length\t"
                 "Mismatches\t"
                 "Gap Openings\t"
                 "Query Start\t"
                 "Query End\t"
                 "Subject Start\t"
                 "Subject Eng\t"
                 "E Value\t"
                 "Coverage\t"
                 "Informativeness\t"
                 "Species\t"
                 "Origin Database\t"
                 "Frame\t"
                 "Seed ortholog\t"
                 "Seed E Value\t"
                 "Seed Score\t"
                 "Tax Score\t"
                 "OGs\t"
                 "GO Terms\t"
                 "KEGG Terms"
         <<std::endl;
    for (auto &pair : SEQUENCES) {
        std::cout<<pair.second.print_eggnog()<<std::endl;
    }
    file.close();
}

