/*
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * 2017
*/

//
//#include "ModInterpro.h"
//#include "../ExceptionHandler.h"
//
//std::pair<bool, std::string> ModInterpro::verify_files() {
//
//    if (_is_overwrite) {
//        boostFS::remove_all(interpro_out_dir);
//    } else {
//        boostFS::path file(_new_input);
//        std::string new_out = interpro_out_dir + file.filename().string() + ".tsv";
//        boostFS::path file2(_input_no_hits);
//        std::string no_hits = interpro_out_dir + file2.filename().string() + ".tsv";
//        if (verify_files(new_out,no_hits)) {
//            out.first = new_out; out.second = no_hits;
//            parse_results_interpro(SEQUENCES,out);
//            return;
//        }
//    }
//    return std::pair<bool, std::string>();
//}
//
//void ModInterpro::execute(std::map<std::string, QuerySequence> &) {
//
//    print_debug("Executing InterProScan...");
//    std::string interpro_out_dir = _outpath + ONTOLOGY_OUT_PATH;
//    std::string annotation_std = interpro_out_dir + "annotation_std";
//    std::pair<std::string,std::string> out;
//
//
//    boostFS::create_directories(interpro_out_dir);
//    std::unordered_map<std::string,std::string> command_map= {
//            {"-i",""},
//            {"-goterms",""},
//            {"-iprlookup",""},
//            {"-pa", ""},
//            {"-d",interpro_out_dir}
//    };
//    int ct = 0;
//    if (!databases.empty()) {
//        command_map["-appl"] = "";
//        for (std::string &val : databases) {
//            if (ct != 0) command_map["-appl"]+=",";
//            command_map["-appl"]+=val;
//            ct++;
//        }
//    }
//    for (int i=0; i<2;i++) {
//        std::string path;
//        i == 0 ? path = _new_input : path = _input_no_hits;
//        if (!file_exists(path)) {
//            print_debug("File not found at: " + path + " skipping...");
//            continue;
//        }
//        std::ifstream inFile(_input_no_hits);
//        long line_num = std::count(std::istreambuf_iterator<char>(inFile),
//                                   std::istreambuf_iterator<char>(), '\n');
//        inFile.close();
//        if (line_num < 2) continue;
//        command_map["-i"] = path;
//        boostFS::path file(_new_input);
//        std::string filename = file.filename().string();
//        std::string cmd = generate_command(command_map,_ontology_exe);
//        print_debug("\nExecuting InterProScan against protein sequences...\n"
//                    + cmd);
//        if (execute_cmd(cmd, annotation_std) !=0) {
//            throw ExceptionHandler("Error executing eggnog mapper", ENTAP_ERR::E_RUN_ANNOTATION);
//        }
//        i == 0 ? out.first=interpro_out_dir + filename+".tsv" : out.second=interpro_out_dir + filename+".tsv";
//    }
//    parse_results_interpro(SEQUENCES,out);
//
//}
//
//void ModInterpro::parse(std::map<std::string, QuerySequence> &) {
//
//
//    std::map<std::string, struct_go_term> GO_DATABASE;
//    std::map<std::string, interpro_struct> interpro_map;
//
//    std::string msg = ENTAP_STATS::SOFTWARE_BREAK + "Ontology - Interpro\n" +
//                      ENTAP_STATS::SOFTWARE_BREAK;
//    print_statistics(msg);
////    const std::string KEY_PROTEIN_DATA      = _HEADERS[0];
////    const std::string KEY_PROTEIN_ID        = _HEADERS[1];
////    const std::string KEY_PROTEIN_TERM      = _HEADERS[2];
////    const std::string KEY_E_VALUE           = _HEADERS[3];
////    const std::string KEY_INTERPRO_ID       = _HEADERS[4];
////    const std::string KEY_INTERPRO_TERM     = _HEADERS[5];
////    const std::string KEY_PATHWAY           = _HEADERS[9];
//
//    try {
//        GO_DATABASE = read_go_map();
//    } catch (ExceptionHandler const &e) {throw e;}
//    for (int i=0; i<2;i++) {
//        std::string path;
//        i == 0 ? path=out.first : path=out.second;
//        print_debug("Interpro file located at " + path + " being filtered");
//        if (!file_exists(path)) {
//            print_debug("File not found, skipping...");continue;
//        }
//        interpro_format_fix(path);
//        std::string qseqid, temp, protein, data_id, data_term, score, score2, temp2,
//                data, ipr_id, ipr_term,go_id,path_id,temp3;
//        double e_val;
//        io::CSVReader<INTERPRO_COL_NUM, io::trim_chars<' '>,
//        io::no_quote_escape<'\t'>> in(out.first);
//        in.set_header("qseqid", "temp", "temp3","protein",
//                      "data_id", "data_term", "score", "score2", "e_val", "temp2",
//                      "data", "ipr_id", "ipr_term","go_id","path_id");
//        while (in.read_row(qseqid, temp, temp3,protein, data_id, data_term, score, score2, e_val, temp2,
//                           data, ipr_id, ipr_term,go_id,path_id)) {
//            std::map<std::string, interpro_struct>::iterator iterator = interpro_map.find(qseqid);
//            if (iterator != interpro_map.end())if (iterator->second._eval < e_val) continue;
//            std::map<std::string,std::string> out_map;
//            interpro_struct out_struct;
////            out_map[KEY_PROTEIN_DATA ]     =    protein;
////            out_map[KEY_PROTEIN_ID   ]     =    data_id;
////            out_map[KEY_PROTEIN_TERM ]     =    data_term;
////            out_map[KEY_E_VALUE      ]     =    std::to_string(e_val);
////            out_map[KEY_INTERPRO_ID  ]     =    ipr_id;
////            out_map[KEY_INTERPRO_TERM]     =    ipr_term;
////            out_map[KEY_PATHWAY      ]     =    path_id;
//
//            out_struct._eval = e_val;
//            out_struct._go_map = parse_go_list(go_id,GO_DATABASE,'|');
//            out_struct._results = out_map;
//            interpro_map[qseqid]=out_struct;
//        }
//    }
//
//    // TODO stats
//    for (auto &pair : SEQUENCES) {
//        std::map<std::string, interpro_struct>::iterator it = interpro_map.find(pair.first);
//        if (it != interpro_map.end()) {
//            pair.second.set_ontology_results(it->second._results);
////            pair.second.set_go_parsed(it->second._go_map);
//        }
//    }
//    GO_DATABASE.clear();
//    print_interpro(SEQUENCES);
//
//}
//
//void ModInterpro::set_data(std::vector<short> &, std::string &, int) {
//
//}
//
//
//
//void Ontology::interpro_format_fix(std::string& path) {
//    std::string out_path = path + "_alt";
//    std::ifstream file(path);
//    std::ofstream out(path+"_alt");
//    std::string line;
//    while (std::getline(file,line)) {
//        if (line.empty()) continue;
//        long ct = INTERPRO_COL_NUM - std::count(line.begin(),line.end(),'\t') - 1;
//        out << line;
//        for (long i = ct ; i >0; i--) out << '\t';
//        out << std::endl;
//    }
//    file.close();out.close();
//    boostFS::remove(path);
//    boostFS::rename(out_path,path);
//}


