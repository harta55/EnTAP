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
#include "SimilaritySearch.h"

namespace boostFS = boost::filesystem;
Ontology::Ontology(int thread, std::string egg_exe, std::string outpath,
                   std::string entap_exe,std::string input,
                   boost::program_options::variables_map &user_input) {
    _ontology_exe = egg_exe;
    _threads = thread;
    _entap_exe = entap_exe;
    _outpath = outpath;
    _new_input = input;
    ONTOLOGY_OUT_PATH = "ontology/";
    _is_overwrite = (bool) user_input.count(ENTAP_CONFIG::INPUT_FLAG_OVERWRITE);
    _software_flag = user_input[ENTAP_CONFIG::INPUT_FLAG_ONTOLOGY].as<short>();
    std::vector<std::string> _interpro_databases =
            user_input[ENTAP_CONFIG::INPUT_FLAG_INTERPRO].as<std::vector<std::string>>();
}


void Ontology::execute(query_map_struct &SEQUENCES, std::string input,std::string no_hit) {
    _new_input = input;
    _input_no_hits = no_hit;
    init_headers();
    try {
        switch(_software_flag) {
            case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
                run_eggnog(SEQUENCES);
                break;
            case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
                run_interpro(SEQUENCES,_interpro_databases);
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
    entapInit:: print_msg("Beginning to parse eggnog results...");
    std::string msg = ENTAP_STATS::SOFTWARE_BREAK + "Ontology - Eggnog\n" +
                      ENTAP_STATS::SOFTWARE_BREAK;
    entapExecute::print_statistics(msg, _outpath);
    DatabaseHelper database;
    if (!database.open(_entap_exe + ENTAP_CONFIG::GO_DB_PATH))
        throw ExceptionHandler("Unable to open GO database",ENTAP_ERR::E_PARSE_EGGNOG);
    std::map<std::string, int> eggnog_map;
    unsigned int count_total_go_hits=0, count_total_go_terms=0, count_no_go=0,count_no_kegg=0,
            count_TOTAL_hits=0, count_total_kegg_terms=0, count_total_kegg_hits=0;
    for (int i=0; i<2;i++) {
        std::string path;
        i == 0 ? path=out.first : path=out.second;
        entapInit::print_msg("Eggnog file located at " + path + " being filtered");
        if (!entapInit::file_exists(path)) {
            entapInit::print_msg("File not found, skipping...");continue;
        }
        path = eggnog_format(path);
        std::string qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                best_og, cog_cat, eggnog_annot;
        io::CSVReader<ENTAP_EXECUTE::EGGNOG_COL_NUM, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(path);
        // io::single_line_comment<'#'>??
        while (in.read_row(qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                           best_og, cog_cat, eggnog_annot)) {
            query_map_struct::iterator it = SEQUENCES.find(qseqid);
            if (it != SEQUENCES.end()) {
                it->second.set_eggnog_results(seed_ortho,seed_e,seed_score,predicted_gene,go_terms,
                                              kegg,tax_scope,ogs);
                it->second.set_go_parsed(parse_go_list(go_terms,database,','));
            }
            count_TOTAL_hits++;
            eggnog_map[qseqid] = 1;
            if (!go_terms.empty()) {
                count_total_go_hits++;
                long ct = std::count(go_terms.begin(), go_terms.end(), ',');
                count_total_go_terms += ct + 1;
            } else {
                count_no_go++;
            }
            if (!kegg.empty()) {
                count_total_kegg_hits++;
                long ct = std::count(kegg.begin(), kegg.end(), ',');
                count_total_kegg_terms += ct + 1;
            } else {
                count_no_kegg++;
            }
        }
        boostFS::remove(path);
    }
    entapInit::print_msg("Success! Computing overall statistics...");
    unsigned int count_no_hits=0;
    for (auto &pair : SEQUENCES) {
        if (eggnog_map.find(pair.first) == eggnog_map.end()) count_no_hits++;
    }
    std::stringstream ss;
    ss <<
       "Statistics for overall Eggnog results: " <<
       "\nTotal hits: " << count_TOTAL_hits <<
       "\nTotal sequences that did not hit against Eggnog databases: " <<count_no_hits<<
       "\nTotal sequences containing GO terms: " << count_total_go_hits<<
       "\nTotal sequences that did not have GO terms: " << count_no_go <<
       "\nTotal matched GO terms: " << count_total_go_terms <<
       "\nTotal sequences containing KEGG terms: " << count_total_kegg_hits<<
       "\nTotal sequences that did not have KEGG terms: " << count_no_kegg<<
       "\nTotal matched KEGG terms: " << count_total_kegg_terms;
    msg = ss.str();
    entapExecute::print_statistics(msg,_outpath);
    database.close();
    entapInit::print_msg("Success!");
    print_eggnog(SEQUENCES);
}


void Ontology::run_eggnog(query_map_struct &SEQUENCES) {
    entapInit::print_msg("Running eggnog...");
    std::string eggnog_out_dir = _outpath + ONTOLOGY_OUT_PATH;
    std::string annotation_base_flag = eggnog_out_dir + "annotation_results";
    std::string annotation_no_flag = eggnog_out_dir + "annotation_results_no_hits";
    std::string annotation_std = eggnog_out_dir + "annotation_std";
    std::string eggnog_command = "python " + _ontology_exe + " ";
    std::pair<std::string,std::string> out;
    if (_is_overwrite) {
        boostFS::remove_all(eggnog_out_dir);
    } else {
        std::string hit_out = annotation_base_flag  +".emapper.annotations";
        std::string no_hit_out = annotation_no_flag +".emapper.annotations";
        if (verify_files(hit_out, no_hit_out)) {
            out.first = hit_out;out.second = no_hit_out;
            parse_results_eggnog(SEQUENCES, out);
            return;
        }
    }
    boostFS::create_directories(eggnog_out_dir);
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
        std::ifstream inFile(_input_no_hits);
        long line_num = std::count(std::istreambuf_iterator<char>(inFile),
                   std::istreambuf_iterator<char>(), '\n');
        inFile.close();
        if (line_num >1) {
            eggnog_command_map["-i"] = _input_no_hits;
            eggnog_command_map["--output"] = annotation_no_flag;
            eggnog_command = "python " + _ontology_exe + " ";
            for (auto &pair : eggnog_command_map) eggnog_command += pair.first + " " + pair.second + " ";
            entapInit::print_msg("\nExecuting eggnog mapper against protein sequences that did not hit databases...\n"
                                 + eggnog_command);
            if (entapInit::execute_cmd(eggnog_command, annotation_std) !=0) {
                throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
            }
            out.second = annotation_no_flag + ".emapper.annotations";
        }
    }
    entapInit::print_msg("Success!");
    parse_results_eggnog(SEQUENCES, out);
}

std::map<std::string,std::vector<std::string>> Ontology::parse_go_list
        (std::string list, DatabaseHelper &database,char delim) {
    std::map<std::string,std::vector<std::string>> output;
    if (list.empty()) return output;
    std::istringstream ss(list);
    std::string temp;
    std::vector<std::vector<std::string>>results;
    while (std::getline(ss,temp,delim)) {
        char *query = sqlite3_mprintf(
                "SELECT category,term from terms WHERE goid=%Q",temp.c_str());
        try {
            results = database.query(query);
            if (!results.empty())output[results[0][0]].push_back(temp + "-" + results[0][1]);
        } catch (std::exception e) {
            throw e;
        }
    }
    return output;
}

void Ontology::print_eggnog(Ontology::query_map_struct &SEQUENCES) {
    entapInit::print_msg("Beginning to print final results...");
    std::string final_annotations = _outpath + "final_annotations.tsv";
    boostFS::remove(final_annotations);
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
                 "Subject End\t"
                 "E Value\t"
                 "Coverage\t"
                 "Informativeness\t"
                 "Species\t"
                 "Origin Database\t"
                 "Frame\t"
                 "Seed ortholog\t"
                 "Seed E Value\t"
                 "Seed Score\t"
                 "Predicted Gene\t"
                 "Tax Scope\t"
                 "OGs\t"
                 "KEGG Terms\t"
                 "GO Biological\t"
                 "GO Cellular\t"
                 "GO Molecular"
         <<std::endl;
    for (auto &pair : SEQUENCES) {
        file<< pair.second.print_final_results(_software_flag,_HEADERS)<<std::endl;
    }
    entapInit::print_msg("Success!");
    file.close();
}

std::string Ontology::eggnog_format(std::string file) {
    std::string out_path = file + "_alt";
    boostFS::remove(out_path);
    std::ifstream in_file(file);
    std::ofstream out_file(out_path);
    std::string line;
    while (getline(in_file,line)) {
        if (line.at(0) == '#' || line.empty()) continue;
        out_file << line << std::endl;
    }
    in_file.close();
    out_file.close();
    return out_path;
}

void Ontology::run_interpro(Ontology::query_map_struct &SEQUENCES, std::vector<std::string>& databases) {
    entapInit::print_msg("Executing InterProScan...");
    std::string interpro_out_dir = _outpath + ONTOLOGY_OUT_PATH;
    std::string annotation_std = interpro_out_dir + "annotation_std";
    std::pair<std::string,std::string> out;

    if (_is_overwrite) {
        boostFS::remove_all(interpro_out_dir);
    } else {
        boostFS::path file(_new_input);
        std::string new_out = interpro_out_dir + file.filename().string() + ".tsv";
        boostFS::path file2(_input_no_hits);
        std::string no_hits = interpro_out_dir + file2.filename().string() + ".tsv";
        if (verify_files(new_out,no_hits)) {
            out.first = new_out; out.second = no_hits;
            parse_results_interpro(SEQUENCES,out);
            return;
        }
    }
    boostFS::create_directories(interpro_out_dir);
    std::unordered_map<std::string,std::string> command_map= {
            {"-i",""},
            {"-goterms",""},
            {"-iprlookup",""},
            {"-pa", ""},
            {"-d",interpro_out_dir}
    };
    int ct = 0;
    if (!databases.empty()) {
        command_map["-appl"] = "";
        for (std::string &val : databases) {
            if (ct != 0) command_map["-appl"]+=",";
            command_map["-appl"]+=val;
            ct++;
        }
    }
    for (int i=0; i<2;i++) {
        std::string path;
        i == 0 ? path = _new_input : path = _input_no_hits;
        if (!entapInit::file_exists(path)) {
            entapInit::print_msg("File not found at: " + path + " skipping...");
            continue;
        }
        std::ifstream inFile(_input_no_hits);
        long line_num = std::count(std::istreambuf_iterator<char>(inFile),
                                   std::istreambuf_iterator<char>(), '\n');
        inFile.close();
        if (line_num < 2) continue;
        command_map["-i"] = path;
        boostFS::path file(_new_input);
        std::string filename = file.filename().string();
        std::string cmd = entapInit::generate_command(command_map,_ontology_exe);
        entapInit::print_msg("\nExecuting InterProScan against protein sequences...\n"
                             + cmd);
        if (entapInit::execute_cmd(cmd, annotation_std) !=0) {
            throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
        }
        i == 0 ? out.first=interpro_out_dir + filename+".tsv" : out.second=interpro_out_dir + filename+".tsv";
    }
    parse_results_interpro(SEQUENCES,out);
}

void Ontology::parse_results_interpro(Ontology::query_map_struct &SEQUENCES,
                                      std::pair<std::string, std::string> &out) {
    std::string msg = ENTAP_STATS::SOFTWARE_BREAK + "Ontology - Interpro\n" +
                      ENTAP_STATS::SOFTWARE_BREAK;
    entapExecute::print_statistics(msg, _outpath);
    DatabaseHelper database;
    const std::string KEY_PROTEIN_DATA      = _HEADERS[0];
    const std::string KEY_PROTEIN_ID        = _HEADERS[1];
    const std::string KEY_PROTEIN_TERM      = _HEADERS[2];
    const std::string KEY_E_VALUE           = _HEADERS[3];
    const std::string KEY_INTERPRO_ID       = _HEADERS[4];
    const std::string KEY_INTERPRO_TERM     = _HEADERS[5];
    const std::string KEY_PATHWAY           = _HEADERS[9];

    if (!database.open(_entap_exe + ENTAP_CONFIG::GO_DB_PATH))
        throw ExceptionHandler("Unable to open GO database",ENTAP_ERR::E_PARSE_EGGNOG);
    struct interpro_struct {
        double _eval;
        std::map<std::string,std::string> _results;
        std::map<std::string,std::vector<std::string>> _go_map;
    };
    std::map<std::string, interpro_struct> interpro_map;
    for (int i=0; i<2;i++) {
        std::string path;
        i == 0 ? path=out.first : path=out.second;
        entapInit::print_msg("Interpro file located at " + path + " being filtered");
        if (!entapInit::file_exists(path)) {
            entapInit::print_msg("File not found, skipping...");continue;
        }
        interpro_format_fix(path);
        std::string qseqid, temp, protein, data_id, data_term, score, score2, temp2,
                data, ipr_id, ipr_term,go_id,path_id,temp3;
        double e_val;
        io::CSVReader<ENTAP_EXECUTE::INTERPRO_COL_NUM, io::trim_chars<' '>,
                io::no_quote_escape<'\t'>> in(out.first);
        in.set_header("qseqid", "temp", "temp3","protein",
                      "data_id", "data_term", "score", "score2", "e_val", "temp2",
                      "data", "ipr_id", "ipr_term","go_id","path_id");
        while (in.read_row(qseqid, temp, temp3,protein, data_id, data_term, score, score2, e_val, temp2,
                           data, ipr_id, ipr_term,go_id,path_id)) {
            std::map<std::string, interpro_struct>::iterator iterator = interpro_map.find(qseqid);
            if (iterator != interpro_map.end())if (iterator->second._eval < e_val) continue;
            std::map<std::string,std::string> out_map;
            interpro_struct out_struct;
            out_map[KEY_PROTEIN_DATA ]     =    protein;
            out_map[KEY_PROTEIN_ID   ]     =    data_id;
            out_map[KEY_PROTEIN_TERM ]     =    data_term;
            out_map[KEY_E_VALUE      ]     =    std::to_string(e_val);
            out_map[KEY_INTERPRO_ID  ]     =    ipr_id;
            out_map[KEY_INTERPRO_TERM]     =    ipr_term;
            out_map[KEY_PATHWAY      ]     =    path_id;

            out_struct._eval = e_val;
            out_struct._go_map = parse_go_list(go_id,database,'|');
            out_struct._results = out_map;
            interpro_map[qseqid]=out_struct;
        }
    }

    // TODO stats
    for (auto &pair : SEQUENCES) {
        std::map<std::string, interpro_struct>::iterator it = interpro_map.find(pair.first);
        if (it != interpro_map.end()) {
            pair.second.set_ontology_results(it->second._results);
            pair.second.set_go_parsed(it->second._go_map);
        }
    }
    database.close();
    print_interpro(SEQUENCES);
}

void Ontology::init_headers() {
    switch (_software_flag) {
        case ENTAP_EXECUTE::EGGNOG_INT_FLAG:
            _HEADERS = {
                    "Seed ortholog"     ,
                    "Seed E Value"      ,
                    "Seed Score"        ,
                    "Predicted Gene"    ,
                    "Tax Scope"         ,
                    "OGs"               ,
                    "GO Biological"     ,
                    "GO Cellular"       ,
                    "GO Molecular"      ,
                    "KEGG Terms"
            };
            break;
        case ENTAP_EXECUTE::INTERPRO_INT_FLAG:
            _HEADERS = {
                "Protein Database"      ,
                "Protein ID"            ,
                "Protein Term"          ,
                "E_value"               ,
                "InterPro ID"           ,
                "InterPro Term"         ,
                "GO Biological"         ,
                "GO Cellular"           ,
                "GO Molecular"          ,
                "Pathway Terms"
            };
            break;
        default:
            break;
    }
}

// TODO combine with eggnog
void Ontology::print_interpro(Ontology::query_map_struct &SEQUENCES) {
    std::string final_annotations = _outpath + "final_annotations.tsv";
    print_header(final_annotations);
    std::ofstream file(final_annotations, std::ios::out | std::ios::app);
    for (auto &pair : SEQUENCES) {
        file<< pair.second.print_final_results(_software_flag,_HEADERS)<<std::endl;
    }
    file.close();
}

void Ontology::print_header(std::string file) {
    std::ofstream ofstream(file, std::ios::out | std::ios::app);
    ofstream <<
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
                     "Frame";
    for (std::string val : _HEADERS) ofstream << '\t' << val;
    ofstream<<std::endl;
    ofstream.close();
}

bool Ontology::verify_files(std::string hits,std::string no_hits) {
    bool verified = false;
    entapInit::print_msg("Overwrite was unselected, verifying output files...");
    if (entapInit::file_exists(hits)) {
        entapInit::print_msg("File located at: " + hits + " found");
        verified = true;
    } else entapInit::print_msg("File located at: " + hits + " NOT found");
    if (entapInit::file_exists(no_hits)) {
            entapInit::print_msg("File located at: " + no_hits + " found");
            verified = true;
    } else entapInit::print_msg("File located at: " + no_hits + " NOT found");
    if (verified) {
        entapInit::print_msg("One or more ontology files were found, skipping ontology execution");
    } else {
        entapInit::print_msg("No ontology files were found, continuing with execution");
    }
    return verified;
}

void Ontology::interpro_format_fix(std::string& path) {
    std::string out_path = path + "_alt";
    std::ifstream file(path);
    std::ofstream out(path+"_alt");
    std::string line;
    while (std::getline(file,line)) {
        if (line.empty()) continue;
        long ct = ENTAP_EXECUTE::INTERPRO_COL_NUM -
                  std::count(line.begin(),line.end(),'\t') - 1;
        out << line;
        for (long i = ct ; i >0; i--) out << '\t';
        out << std::endl;
    }
    file.close();out.close();
    boostFS::remove(path);
    boostFS::rename(out_path,path);
}
