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

namespace boostFS = boost::filesystem;
Ontology::Ontology(int thread, std::string egg_exe, std::string outpath, std::string entap_exe,
                   std::string input, boost::program_options::variables_map &user_input,
                   std::string database, GraphingManager* graphing) {
    print_debug("Spawn object - Ontology");
    _ontology_exe = egg_exe;
    _threads = thread;
    _entap_exe = entap_exe;
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
    _eggnog_db_path = database;
    _graphingManager = graphing;
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
    print_debug("Beginning to parse eggnog results...");

    typedef std::map<std::string,std::map<std::string, unsigned int>> GO_top_map_t;

    std::stringstream                        ss;
    std::string                              out_msg;
    std::string                              out_no_hits_nucl;
    std::string                              out_no_hits_prot;
    std::string                              out_hit_nucl;
    std::string                              out_hit_prot;
    std::string                              path;
    std::string                              fig_txt_bar_go_overall;
    std::string                              fig_png_bar_go_overall;
    std::string                              fig_txt_bar_ortho;
    std::string                              fig_png_bar_ortho;
    std::string                              tax_scope_readable;
    std::map<std::string, struct_go_term>    GO_DATABASE;
    std::map<std::string, int>               eggnog_map;
    unsigned int                             count_total_go_hits=0;
    unsigned int                             count_total_go_terms=0;
    unsigned int                             count_go_bio=0;
    unsigned int                             count_go_cell=0;
    unsigned int                             count_go_mole=0;
    unsigned int                             count_no_go=0;
    unsigned int                             count_no_kegg=0;
    unsigned int                             count_TOTAL_hits=0;         // All ortho matches
    unsigned int                             count_total_kegg_terms=0;
    unsigned int                             count_total_kegg_hits=0;
    unsigned int                             count_no_hits=0;            // Unannotated OGs
    unsigned int                             count_tax_scope=0;
    unsigned int                             ct = 0;
    float                                    percent;
    DatabaseHelper                           EGGNOG_DATABASE;
    std::map<std::string, unsigned int>      tax_scope_ct_map;
    GO_top_map_t                             go_combined_map;     // Just for convenience
    go_struct                                go_parsed;
    GraphingStruct                           graphingStruct;

    boostFS::remove_all(_processed_dir);
    boostFS::create_directories(_processed_dir);

    ss<<std::fixed<<std::setprecision(2);

    try {
        GO_DATABASE = read_go_map();
    } catch (ExceptionHandler const &e) {throw e;}

    if (!EGGNOG_DATABASE.open(_eggnog_db_path))
        throw ExceptionHandler("Unable to open GO database",ENTAP_ERR::E_PARSE_EGGNOG);

    for (int i=0; i<2;i++) {
        i == 0 ? path=out.first : path=out.second;
        print_debug("Eggnog file located at " + path + " being filtered");
        if (!file_exists(path)) {
            print_debug("File not found, skipping...");continue;
        }
        path = eggnog_format(path);
        std::string qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                best_og, cog_cat, eggnog_annot;
        io::CSVReader<EGGNOG_COL_NUM, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(path);
        while (in.read_row(qseqid, seed_ortho, seed_e, seed_score, predicted_gene, go_terms, kegg, tax_scope, ogs,
                           best_og, cog_cat, eggnog_annot)) {
            query_map_struct::iterator it = SEQUENCES.find(qseqid);
            if (it != SEQUENCES.end()) {
                count_TOTAL_hits++;
                it->second.set_eggnog_results(seed_ortho,seed_e,seed_score,predicted_gene,go_terms,
                                              kegg,tax_scope,ogs, EGGNOG_DATABASE);
                go_parsed = parse_go_list(go_terms,GO_DATABASE,',');
                it->second.set_go_parsed(go_parsed, ENTAP_EXECUTE::EGGNOG_INT_FLAG);
                it->second.set_is_family_assigned(true);
                eggnog_map[qseqid] = 1;
                if (!go_parsed.empty()) {
                    count_total_go_hits++;
                    it->second.set_is_one_go(true);
                    for (auto &pair : go_parsed) {
                        for (std::string &term : pair.second) {
                            count_total_go_terms++;
                            if (pair.first == GO_MOLECULAR_FLAG) {
                                count_go_mole++;
                            } else if (pair.first == GO_CELLULAR_FLAG) {
                                count_go_cell++;
                            } else if (pair.first == GO_BIOLOGICAL_FLAG) {
                                count_go_bio++;
                            }
                            if (go_combined_map[pair.first].count(term)) {
                                go_combined_map[pair.first][term]++;
                            } else go_combined_map[pair.first][term] = 1;
                            if (go_combined_map[GO_OVERALL_FLAG].count(term)) {
                                go_combined_map[GO_OVERALL_FLAG][term]++;
                            } else go_combined_map[GO_OVERALL_FLAG][term]=1;
                        }
                    }
                } else {
                    count_no_go++;
                }

                // Compile KEGG information
                if (!kegg.empty()) {
                    count_total_kegg_hits++;
                    ct = (unsigned int) std::count(kegg.begin(), kegg.end(), ',');
                    count_total_kegg_terms += ct + 1;
                    it->second.set_is_one_kegg(true);
                } else {
                    count_no_kegg++;
                }

                // Compile Taxonomic Orthogroup stats
                tax_scope_readable = it->second.get_tax_scope();
                if (!tax_scope_readable.empty()) {
                    count_tax_scope++;
                    if (tax_scope_ct_map.count(tax_scope_readable)) {
                        tax_scope_ct_map[tax_scope_readable]++;
                    } else tax_scope_ct_map[tax_scope_readable] = 1;

                }
            }
        }
        boostFS::remove(path);
    }

    EGGNOG_DATABASE.close();
    out_no_hits_nucl = (boostFS::path(_processed_dir) / boostFS::path(OUT_UNANNOTATED_NUCL)).string();
    out_no_hits_prot = (boostFS::path(_processed_dir) / boostFS::path(OUT_UNANNOTATED_PROT)).string();
    out_hit_nucl     = (boostFS::path(_processed_dir) / boostFS::path(OUT_ANNOTATED_NUCL)).string();
    out_hit_prot     = (boostFS::path(_processed_dir) / boostFS::path(OUT_ANNOTATED_PROT)).string();
    std::ofstream file_no_hits_nucl(out_no_hits_nucl, std::ios::out | std::ios::app);
    std::ofstream file_no_hits_prot(out_no_hits_prot, std::ios::out | std::ios::app);
    std::ofstream file_hits_nucl(out_hit_nucl, std::ios::out | std::ios::app);
    std::ofstream file_hits_prot(out_hit_prot, std::ios::out | std::ios::app);

    print_debug("Success! Computing overall statistics...");
    for (auto &pair : SEQUENCES) {
        if (eggnog_map.find(pair.first) == eggnog_map.end()) {
            // Unannotated sequence
            if (!pair.second.get_sequence_n().empty()) file_no_hits_nucl<<pair.second.get_sequence_n()<<std::endl;
            file_no_hits_prot << pair.second.get_sequence_p() << std::endl;
            count_no_hits++;
        } else {
            // Annotated sequence
            if (!pair.second.get_sequence_n().empty()) file_hits_nucl<<pair.second.get_sequence_n()<<std::endl;
            file_hits_prot << pair.second.get_sequence_p() << std::endl;
        }
    }

    file_hits_nucl.close();
    file_hits_prot.close();
    file_no_hits_nucl.close();
    file_no_hits_prot.close();

    ss << ENTAP_STATS::SOFTWARE_BREAK + "Ontology - Eggnog\n" +
          ENTAP_STATS::SOFTWARE_BREAK            <<
       "Statistics for overall Eggnog results: " <<
       "\nTotal sequences with family assignment: " << count_TOTAL_hits <<
       "\nTotal sequences without family assignment: " <<count_no_hits;

    // -------- Top Ten Taxonomic Scopes ------- //
    if (!tax_scope_ct_map.empty()) {
        std::string fig_txt_tax_bar = (boostFS::path(_figure_dir) / GRAPH_EGG_TAX_BAR_TXT).string();
        std::string fig_png_tax_bar = (boostFS::path(_figure_dir) / GRAPH_EGG_TAX_BAR_PNG).string();
        std::ofstream file_tax_bar(fig_txt_tax_bar, std::ios::out | std::ios::app);
        file_tax_bar << "Taxonomic Scope\tCount" << std::endl;

        ss << "\nTop 10 Taxonomic Scopes Assigned:";
        ct = 1;
        std::vector<count_pair> tax_scope_vect(tax_scope_ct_map.begin(),tax_scope_ct_map.end());
        std::sort(tax_scope_vect.begin(),tax_scope_vect.end(),compair());
        for (count_pair &pair : tax_scope_vect) {
            if (ct > 10) break;
            percent = ((float)pair.second / count_tax_scope) * 100;
            ss <<
               "\n\t" << ct << ")" << pair.first << ": " << pair.second <<
               "(" << percent << "%)";
            file_tax_bar << pair.first << '\t' << std::to_string(pair.second) << std::endl;
            ct++;
        }
        file_tax_bar.close();
        graphingStruct.fig_out_path = fig_png_tax_bar;
        graphingStruct.text_file_path = fig_txt_tax_bar;
        graphingStruct.graph_title = GRAPH_EGG_TAX_BAR_TITLE;
        graphingStruct.software_flag = GRAPH_ONTOLOGY_FLAG;
        graphingStruct.graph_type = GRAPH_TOP_BAR_FLAG;
        _graphingManager->graph(graphingStruct);
    }
    // --------------------------------------- //

    ss<<
      "\nTotal sequences with at least one GO term: " << count_total_go_hits <<
      "\nTotal sequences without GO terms: " << count_no_go <<
      "\nTotal GO terms assigned: " << count_total_go_terms;

    if (count_total_go_hits > 0) {
        std::map<std::string, unsigned int> overall_go_tops;
        for (auto &pair : go_combined_map) {
            // Count maps (biological/molecular/cellular/overall)
            std::string fig_txt_go_bar = (boostFS::path(_figure_dir) / pair.first).string() + GRAPH_GO_END_TXT;
            std::string fig_png_go_bar = (boostFS::path(_figure_dir) / pair.first).string() + GRAPH_GO_END_PNG;
            std::ofstream file_go_bar(fig_txt_go_bar, std::ios::out | std::ios::app);
            std::vector<count_pair> go_vect(pair.second.begin(),pair.second.end());
            std::sort(go_vect.begin(),go_vect.end(),compair());
            file_go_bar << "Gene Ontology Term\tCount" << std::endl;
            ss << "\nTop 10" << pair.first << " terms assigned: ";
            ct = 1;
            for (count_pair &pair2 : go_vect) {
                if (ct > 10) break;
                percent = ((float)pair2.second / count_total_go_terms) * 100;
                ss <<
                   "\n\t" << ct << ")" << pair2.first << ": " << pair2.second <<
                   "(" << percent << "%)";
                file_go_bar << pair2.first << '\t' << std::to_string(pair2.second) << std::endl;
                ct++;
            }
            file_go_bar.close();
            graphingStruct.fig_out_path = fig_png_go_bar;
            graphingStruct.text_file_path = fig_txt_go_bar;
            if (pair.first == GO_BIOLOGICAL_FLAG) graphingStruct.graph_title = GRAPH_GO_BAR_BIO_TITLE;
            if (pair.first == GO_CELLULAR_FLAG) graphingStruct.graph_title = GRAPH_GO_BAR_CELL_TITLE;
            if (pair.first == GO_MOLECULAR_FLAG) graphingStruct.graph_title = GRAPH_GO_BAR_MOLE_TITLE;
            if (pair.first == GO_OVERALL_FLAG) graphingStruct.graph_title = GRAPH_GO_BAR_ALL_TITLE;
            // Other params can stay the same
            _graphingManager->graph(graphingStruct);
        }
    }
    ss<<
      "\nTotal sequences with at least one pathway (KEGG) assignment: " << count_total_kegg_hits<<
      "\nTotal sequences without pathways (KEGG): " << count_no_kegg<<
      "\nTotal pathways (KEGG) assigned: " << count_total_kegg_terms;
    out_msg = ss.str();
    print_statistics(out_msg);
    GO_DATABASE.clear();
    print_debug("Success!");
    print_eggnog(SEQUENCES);
}


void Ontology::run_eggnog(query_map_struct &SEQUENCES) {
    print_debug("Running eggnog...");

    std::string                        annotation_base_flag;
    std::string                        annotation_no_flag;
    std::string                        annotation_std;
    std::string                        eggnog_command;
    std::string                        hit_out;
    std::string                        no_hit_out;
    std::pair<std::string,std::string> out;


    annotation_base_flag = (boostFS::path(_ontology_dir) / boostFS::path("annotation_results")).string();
    annotation_no_flag   = (boostFS::path(_ontology_dir) / boostFS::path("annotation_results_no_hits")).string();
    annotation_std       = (boostFS::path(_ontology_dir) / boostFS::path("annotation_std")).string();
    eggnog_command       = "python " + _ontology_exe + " ";

    if (_is_overwrite) {
        boostFS::remove_all(_ontology_dir);
    } else {
        hit_out = annotation_base_flag  +".emapper.annotations";
        no_hit_out = annotation_no_flag +".emapper.annotations";
        if (verify_files(hit_out, no_hit_out)) {
            out.first = hit_out;out.second = no_hit_out;
            parse_results_eggnog(SEQUENCES, out);
            return;
        }
    }
    boostFS::create_directories(_ontology_dir);
    std::unordered_map<std::string,std::string> eggnog_command_map = {
            {"-i",_new_input},
            {"--output",annotation_base_flag},
            {"--cpu",std::to_string(_threads)},
            {"-m", "diamond"}
    };
    if (file_exists(_new_input)) {
        for (auto &pair : eggnog_command_map)eggnog_command += pair.first + " " + pair.second + " ";
        print_debug("\nExecuting eggnog mapper against protein sequences that hit databases...\n"
                             + eggnog_command);
        if (execute_cmd(eggnog_command, annotation_std) !=0) {
            throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
        }
        print_debug("Success! Results written to: " + annotation_base_flag);
        out.first = annotation_base_flag + ".emapper.annotations";
    } else {
        throw ExceptionHandler("No input file found at: " + _new_input,
                               ENTAP_ERR::E_RUN_EGGNOG);
    }
    if (file_exists(_input_no_hits)) {
        std::ifstream inFile(_input_no_hits);
        long line_num = std::count(std::istreambuf_iterator<char>(inFile),
                   std::istreambuf_iterator<char>(), '\n');
        inFile.close();
        if (line_num >1) {
            eggnog_command_map["-i"] = _input_no_hits;
            eggnog_command_map["--output"] = annotation_no_flag;
            eggnog_command = "python " + _ontology_exe + " ";
            for (auto &pair : eggnog_command_map) eggnog_command += pair.first + " " + pair.second + " ";
            print_debug("\nExecuting eggnog mapper against protein sequences that did not hit databases...\n"
                                 + eggnog_command);
            if (execute_cmd(eggnog_command, annotation_std) !=0) {
                throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
            }
            out.second = annotation_no_flag + ".emapper.annotations";
        }
    }
    print_debug("Success!");
    parse_results_eggnog(SEQUENCES, out);
}

std::map<std::string,std::vector<std::string>> Ontology::parse_go_list
        (std::string list, std::map<std::string,struct_go_term> &GO_DATABASE,char delim) {

    std::map<std::string,std::vector<std::string>> output;
    std::string temp;
    std::vector<std::vector<std::string>>results;

    if (list.empty()) return output;
    std::istringstream ss(list);
    while (std::getline(ss,temp,delim)) {
        struct_go_term term_info = GO_DATABASE[temp];
        output[term_info.category].push_back(temp + "-" + term_info.term +
            "(L=" + term_info.level + ")");
    }
    return output;
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
            *file_map[lvl] << &header << '\t';
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


// TODO remove
std::string Ontology::eggnog_format(std::string file) {

    std::string out_path;
    std::string line;

    out_path = file + "_alt";
    boostFS::remove(out_path);
    std::ifstream in_file(file);
    std::ofstream out_file(out_path);
    while (getline(in_file,line)) {
        if (line.at(0) == '#' || line.empty()) continue;
        out_file << line << std::endl;
    }
    in_file.close();
    out_file.close();
    return out_path;
}

void Ontology::run_interpro(query_map_struct &SEQUENCES, std::vector<std::string>& databases) {
    print_debug("Executing InterProScan...");
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
        if (!file_exists(path)) {
            print_debug("File not found at: " + path + " skipping...");
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
        std::string cmd = generate_command(command_map,_ontology_exe);
        print_debug("\nExecuting InterProScan against protein sequences...\n"
                             + cmd);
        if (execute_cmd(cmd, annotation_std) !=0) {
            throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
        }
        i == 0 ? out.first=interpro_out_dir + filename+".tsv" : out.second=interpro_out_dir + filename+".tsv";
    }
    parse_results_interpro(SEQUENCES,out);
}

void Ontology::parse_results_interpro(query_map_struct &SEQUENCES,
                                      std::pair<std::string, std::string> &out) {

    std::map<std::string, struct_go_term> GO_DATABASE;
    std::map<std::string, interpro_struct> interpro_map;

    std::string msg = ENTAP_STATS::SOFTWARE_BREAK + "Ontology - Interpro\n" +
                      ENTAP_STATS::SOFTWARE_BREAK;
    print_statistics(msg);
//    const std::string KEY_PROTEIN_DATA      = _HEADERS[0];
//    const std::string KEY_PROTEIN_ID        = _HEADERS[1];
//    const std::string KEY_PROTEIN_TERM      = _HEADERS[2];
//    const std::string KEY_E_VALUE           = _HEADERS[3];
//    const std::string KEY_INTERPRO_ID       = _HEADERS[4];
//    const std::string KEY_INTERPRO_TERM     = _HEADERS[5];
//    const std::string KEY_PATHWAY           = _HEADERS[9];

    try {
        GO_DATABASE = read_go_map();
    } catch (ExceptionHandler const &e) {throw e;}
    for (int i=0; i<2;i++) {
        std::string path;
        i == 0 ? path=out.first : path=out.second;
        print_debug("Interpro file located at " + path + " being filtered");
        if (!file_exists(path)) {
            print_debug("File not found, skipping...");continue;
        }
        interpro_format_fix(path);
        std::string qseqid, temp, protein, data_id, data_term, score, score2, temp2,
                data, ipr_id, ipr_term,go_id,path_id,temp3;
        double e_val;
        io::CSVReader<INTERPRO_COL_NUM, io::trim_chars<' '>,
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
//            out_map[KEY_PROTEIN_DATA ]     =    protein;
//            out_map[KEY_PROTEIN_ID   ]     =    data_id;
//            out_map[KEY_PROTEIN_TERM ]     =    data_term;
//            out_map[KEY_E_VALUE      ]     =    std::to_string(e_val);
//            out_map[KEY_INTERPRO_ID  ]     =    ipr_id;
//            out_map[KEY_INTERPRO_TERM]     =    ipr_term;
//            out_map[KEY_PATHWAY      ]     =    path_id;

            out_struct._eval = e_val;
            out_struct._go_map = parse_go_list(go_id,GO_DATABASE,'|');
            out_struct._results = out_map;
            interpro_map[qseqid]=out_struct;
        }
    }

    // TODO stats
    for (auto &pair : SEQUENCES) {
        std::map<std::string, interpro_struct>::iterator it = interpro_map.find(pair.first);
        if (it != interpro_map.end()) {
            pair.second.set_ontology_results(it->second._results);
//            pair.second.set_go_parsed(it->second._go_map);
        }
    }
    GO_DATABASE.clear();
    print_interpro(SEQUENCES);
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

bool Ontology::verify_files(std::string hits,std::string no_hits) {
    bool verified = false;
    print_debug("Overwrite was unselected, verifying output files...");
    if (file_exists(hits)) {
        print_debug("File located at: " + hits + " found");
        verified = true;
    } else print_debug("File located at: " + hits + " NOT found");
    if (file_exists(no_hits)) {
            print_debug("File located at: " + no_hits + " found");
            verified = true;
    } else print_debug("File located at: " + no_hits + " NOT found");
    if (verified) {
        print_debug("One or more ontology files were found, skipping ontology execution");
    } else {
        print_debug("No ontology files were found, continuing with execution");
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
        long ct = INTERPRO_COL_NUM - std::count(line.begin(),line.end(),'\t') - 1;
        out << line;
        for (long i = ct ; i >0; i--) out << '\t';
        out << std::endl;
    }
    file.close();out.close();
    boostFS::remove(path);
    boostFS::rename(out_path,path);
}

std::map<std::string,struct_go_term> Ontology::read_go_map () {
    std::map<std::string,struct_go_term> new_map;
    std::string go_db_path = _entap_exe + ENTAP_CONFIG::GO_DB_PATH;
    try {
        {
            std::ifstream ifs(go_db_path);
            boost::archive::binary_iarchive ia(ifs);
            ia >> new_map;
        }
    } catch (std::exception &exception) {
        throw ExceptionHandler(exception.what(), ENTAP_ERR::E_INIT_GO_SETUP);
    }
    return new_map;
};
