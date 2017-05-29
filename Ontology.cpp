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
#include "DatabaseHelper.h"

namespace boostFS = boost::filesystem;

Ontology::Ontology(int thread, bool overwrite, std::string egg_exe, std::string outpath,
                   std::string entap_exe,std::string input) {
    _eggnog_exe = egg_exe;
    _threads = thread;
    _is_overwrite = overwrite;
    _entap_exe = entap_exe;
    _outpath = outpath;
    _new_input = input;
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
    std::string msg = ENTAP_STATS::SOFTWARE_BREAK + "Ontology - Eggnog\n" +
                      ENTAP_STATS::SOFTWARE_BREAK;
    entapExecute::print_statistics(msg, _outpath);
    DatabaseHelper database;
    if (!database.open(_entap_exe + ENTAP_CONFIG::GO_DB_PATH))
        throw ExceptionHandler("Unable to open GO database",ENTAP_ERR::E_PARSE_EGGNOG);

    for (int i=0; i<2;i++) {
        std::string path;
        i == 0 ? path=out.first : path=out.second;
        entapInit::print_msg("Eggnog file located at " + path + " being filtered");
        if (!entapInit::file_exists(path)) {
            entapInit::print_msg("File not found, skipping...");continue;
        }
        std::string qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                best_og, cog_cat, eggnog_annot;
        unsigned int count_total_go_hits=0, count_total_go_terms=0, count_no_go=0,count_no_kegg=0,
            count_TOTAL_hits=0, count_total_kegg_terms=0, count_total_kegg_hits=0;
        io::CSVReader<ENTAP_EXECUTE::EGGNOG_COL_NUM, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(out.first);
        std::map<std::string, int> eggnog_map;
        in.next_line();in.next_line();in.next_line();
        while (in.read_row(qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                           best_og, cog_cat, eggnog_annot)) {
            SEQUENCES[qseqid].set_eggnog_results(seed_ortho,seed_e,seed_score,predicted_gene,go_terms,
                                                 kegg,tax_scope,ogs);
            count_TOTAL_hits++;
            eggnog_map[qseqid] = 1;
            SEQUENCES[qseqid].set_go_parsed(parse_go_list(go_terms,database));
            if (!go_terms.empty()) {
                count_total_go_hits++;
                long ct = std::count(go_terms.begin(), go_terms.end(), ',');
                count_total_go_terms += ct + 1;
            } else {
                count_no_go++;
            }
            if (!kegg.empty()) {
                count_total_kegg_hits++;
                long ct = std::count(go_terms.begin(), go_terms.end(), ',');
                count_total_kegg_terms += ct + 1;
            } else {
                count_no_kegg++;
            }
        }
        unsigned int count_no_hits=0;
        for (auto &pair : SEQUENCES) {
            if (eggnog_map.find(pair.first) == eggnog_map.end()) count_no_hits++;
        }
        std::stringstream ss;
        ss <<
            "Statistics for file located at: " << path <<
            "\nTotal hits: " << count_TOTAL_hits <<
            "\nTotal sequences that did not hit Eggnog Databases: "<<count_no_hits<<
            "\nTotal sequences containing GO terms: " << count_total_go_hits<<
            "\nTotal sequences that did not have GO terms: " << count_no_go <<
            "\nTotal matched GO terms: " << count_total_go_terms <<
            "\nTotal sequences containing KEGG terms: " << count_total_kegg_hits<<
            "\nTotal sequences that did not have KEGG terms: " << count_no_kegg<<
            "\nTotal matched KEGG terms: " << count_total_kegg_terms;
        msg = ss.str();
        entapExecute::print_statistics(msg,_outpath);
    }
    database.close();
    print_eggnog(SEQUENCES);
}


void Ontology::run_eggnog(query_map_struct &SEQUENCES) {

    std::string eggnog_out_dir = _outpath + ENTAP_EXECUTE::ONTOLOGY_OUT_PATH;
    boostFS::remove_all(eggnog_out_dir);
    boostFS::create_directories(eggnog_out_dir);
//    std::string annotation_base_flag = eggnog_out_dir + "annotation_results";
//    std::string annotation_no_flag = eggnog_out_dir + "annotation_results_no_hits";
        std::string annotation_base_flag = "annotation_results";
    std::string annotation_no_flag = "annotation_results_no_hits";
    std::string annotation_std = eggnog_out_dir + "annotation_std";
    std::string eggnog_command = "python " + _eggnog_exe + " ";
    std::pair<std::string,std::string> out;

    // TODO eggnog overwrite

    std::unordered_map<std::string,std::string> eggnog_command_map = {
            {"-i",_new_input},
            {"--output",annotation_base_flag},
            {"--cpu",std::to_string(_threads)},
            {"-m", "diamond"}
    };
    if (entapInit::file_exists(_new_input)) {
        for (auto &pair : eggnog_command_map)eggnog_command += pair.first + " " + pair.second + " ";
        entapInit::print_msg("\nExecuting eggnog mapper against protein sequences that hit databases...\n"
                             + eggnog_command);
        if (entapInit::execute_cmd(eggnog_command, annotation_std) !=0) {
            throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
        }
        entapInit::print_msg("Success! Results written to: " + annotation_base_flag);
        out.first = annotation_base_flag + ".emapper.annotations";
    } else {
        throw ExceptionHandler("No input file found at: " + _new_input,
                               ENTAP_ERR::E_RUN_EGGNOG);
    }
    if (entapInit::file_exists(_input_no_hits)) {
        std::ifstream inFile("file");
        long line_num = std::count(std::istreambuf_iterator<char>(inFile),
                   std::istreambuf_iterator<char>(), '\n');
        if (line_num >1) {
            eggnog_command_map["-i"] = _input_no_hits;
            eggnog_command_map["--output"] = annotation_no_flag;
            eggnog_command = "python " + _eggnog_exe + " ";
            for (auto &pair : eggnog_command_map) eggnog_command += pair.first + " " + pair.second + " ";
            entapInit::print_msg("\nExecuting eggnog mapper against protein sequences that did not hit databases...\n"
                                 + eggnog_command);
            if (entapInit::execute_cmd(eggnog_command, annotation_std) !=0) {
                throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
            }
            out.second = annotation_no_flag + ".emapper.annotations";
        }
    }
    parse_results_eggnog(SEQUENCES, out);
}

std::map<std::string,std::vector<std::string>> Ontology::parse_go_list
        (std::string list, DatabaseHelper &database) {
    std::map<std::string,std::vector<std::string>> output {

    };
    if (list.empty()) return output;
    std::stringstream ss(list);
    std::string temp;
    std::vector<std::vector<std::string>>results;
    while (ss >> temp) {
        char *query = sqlite3_mprintf(
                "SELECT category,term from terms WHERE goid=%Q",temp.c_str());
        try {
            results = database.query(query);
            output[results[0][0]].push_back(temp + "-" + results[0][1]);
        } catch (ExceptionHandler &e) {throw e;}
        if (ss.peek() == ',')
            ss.ignore();
    }
    return output;
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
                 "GO Biological\t"
                 "GO Cellular\t"
                 "GO Molecular\t"
                 "KEGG Terms"
         <<std::endl;
    for (auto &pair : SEQUENCES) {
        std::cout<<pair.second.print_eggnog()<<std::endl;
    }
    file.close();
}

