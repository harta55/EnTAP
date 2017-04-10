//
// Created by harta on 3/4/17.
//

#include <boost/serialization/unordered_map.hpp>
#include <iostream>
#include <boost/archive/binary_iarchive.hpp>
#include <fstream>
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapConsts.h"
#include "EntapInit.h"
#include <thread>
#include <list>
#include "csv.h"
#include "QuerySequence.h"
#include <boost/regex.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options/variables_map.hpp>

namespace boostFS = boost::filesystem;

namespace entapExecute {

    enum ExecuteStates {

        INIT                = 0x01,
        FRAME_SELECTION     = 0x02,
        RSEM                = 0x04,
        FILTER              = 0x08,
        DIAMOND_RUN         = 0x16,
        DIAMOND_PARSE       = 0x32,
        EXIT                = 0x64
    };

    void execute_main(boost::program_options::variables_map &user_input) {
        entapInit::print_msg("enTAP Executing...");
        boostFS::path working_dir(boostFS::current_path());

        unsigned int supported_threads = std::thread::hardware_concurrency();
        int threads;
        if (user_input["threads"].as<int>() > supported_threads) {
            entapInit::print_msg("Specified thread number is larger than available threads,"
                                         "setting threads to " + std::to_string(supported_threads));
            threads = supported_threads;
        } else {
            threads = user_input["threads"].as<int>();
        }

        std::string current_state = user_input["state"].as<std::string>();
        ExecuteStates state = INIT;
        std::vector<std::string> other_databases;
        if (user_input.count("database")) {
            other_databases = user_input["database"].as<std::vector<std::string>>();
        } else other_databases.push_back(ENTAP_CONFIG::NCBI_NULL);
        std::list<std::string> databases = verify_databases(user_input["uniprot"].as<std::vector<std::string>>(),
            user_input["ncbi"].as<std::vector<std::string>>(),other_databases);

        std::string input_path, rsem_out, genemark_out;
        input_path = user_input["input"].as<std::string>();        // Gradually changes between runs
        while (state != EXIT) {
            try {
                switch (state) {
                    case FRAME_SELECTION:
                        genemark_out = genemarkST(input_path);
                        break;
                    case RSEM:
//                bool is_paired = (bool)user_input.count("paired-end");
//                rsem_out = rsem(input_path,user_input["align"].as<std::string>(),is_paired,threads);
                        break;
                    case FILTER:
//                input_path = filter_transcriptome(genemark_out,rsem_out,user_input["fpkm"].as<float>(),input_path);
                        break;
                    case DIAMOND_RUN:
                        diamond_run(databases,input_path,threads);
                        break;
                    default:
                        break;
                }
                verify_state(current_state);
                state = EXIT;

            } catch (ExceptionHandler &e) {
                throw e;
            }
        }
    }

    std::list<std::string> verify_databases(std::vector<std::string> uniprot, std::vector<std::string>ncbi,
                                            std::vector<std::string> database) {
        entapInit::print_msg("Verifying databases...");
        // return file paths
        std::list<std::string> file_paths;

        std::string path;
        entapInit::print_msg("Verifying uniprot databases...");
        if (uniprot.size()>0) {
            for (auto const& u_flag:uniprot) {
                if (u_flag.compare(ENTAP_CONFIG::INPUT_UNIPROT_NULL)!=0) {
                    path = ENTAP_CONFIG::UNIPROT_INDEX_PATH + u_flag + ".dmnd";
                    if (!entapInit::file_exists(path))
                        throw ExceptionHandler("Database located at: "+path+" not found", ENTAP_ERR::E_INPUT_PARSE);
                    file_paths.push_back(path);
                }
            }
        }
        entapInit::print_msg("Complete");
        entapInit::print_msg("Verifying NCBI databases...");

        if (ncbi.size()>0) {
            for (auto const& u_flag:ncbi) {
                if (u_flag.compare(ENTAP_CONFIG::NCBI_NULL)!=0) {
                    path = ENTAP_CONFIG::NCBI_INDEX_PATH + u_flag + ".dmnd";
                    if (!entapInit::file_exists(path))
                        throw ExceptionHandler("Database located at: "+path+" not found", ENTAP_ERR::E_INPUT_PARSE);
                    file_paths.push_back(path);
                }
            }
        }
        entapInit::print_msg("Complete");
        entapInit::print_msg("Verifying other databases...");
        if (database.size()>0) {
            for (auto const& data_path:database) {
                if (data_path.compare(ENTAP_CONFIG::NCBI_NULL)==0) continue;
                if (!entapInit::file_exists(data_path)) {
                    throw ExceptionHandler("Database located at: "+path+" not found", ENTAP_ERR::E_INPUT_PARSE);
                }
                boostFS::path bpath(data_path); std::string ext = bpath.extension().string();
                if (ext.compare(".dmnd") == 0) {
                    entapInit::print_msg("User has input a diamond indexed database at: " +
                        data_path);
                    file_paths.push_back(data_path); continue;
                } else {
                    entapInit::print_msg("User has input a database at: "+data_path);
                    std::string test_path = ENTAP_CONFIG::BIN_PATH + data_path + ".dmnd";
                    entapInit::print_msg("Checking if indexed file exists at: " + test_path);
                    if (!entapInit::file_exists(test_path)) {
                        throw ExceptionHandler("Database located at: "+path+" not found", ENTAP_ERR::E_INPUT_PARSE);
                    } else {
                        file_paths.push_back(test_path);
                    }
                }
            }
        }
        return file_paths;
    }

    std::string genemarkST(std::string file_path) {
        // Outfiles: file/path.faa, file/path.fnn
        // assumes working directory right now
        entapInit::print_msg("Running genemark...");
        std::string genemark_cmd = ENTAP_EXECUTE::GENEMARK_EXE_PATH + " -faa -fnn " + file_path;
        if (entapInit::execute_cmd(genemark_cmd) != 0 ) {
            throw ExceptionHandler("Error in running genemark at file located at: " +
                file_path, ENTAP_ERR::E_INIT_INDX_DATA_NOT_FOUND);
        }
        entapInit::print_msg("Success!");
        // Format genemarks-t output (remove blank lines)
        entapInit::print_msg("Formatting genemark files");

        std::string genemark_out_dir = ENTAP_EXECUTE::EXECUTE_OUT_PATH + "genemark/";
        boostFS::create_directories(genemark_out_dir);
        boost::filesystem::path file_name(file_path); file_name = file_name.filename();
        std::list<std::string> out_names {file_name.string()+".faa", file_name.string()+".fnn"};
        std::string line;
        for (std::string path : out_names) {
            std::ifstream in_file(path);
            std::string temp_name = path+"_alt";
            std::string out_path = genemark_out_dir+path;
            std::ofstream out_file(path+"_alt");
            while (getline(in_file,line)){
                if (!line.empty()) {
                    out_file << line << '\n';
                }
            }
            in_file.close();
            out_file.close();
            if (remove(path.c_str())!=0 || rename(temp_name.c_str(),out_path.c_str())!=0) {
                throw ExceptionHandler("Error formatting/moving genemark results", ENTAP_ERR::E_INIT_TAX_READ);
            }
        }
        std::string lst_file = file_name.string() + ".lst";
        std::string out_lst = genemark_out_dir + lst_file;

        if (rename(lst_file.c_str(),out_lst.c_str())!=0) {
            throw ExceptionHandler("Error moving genemark results", ENTAP_ERR::E_INIT_TAX_READ);
        }
        return genemark_out_dir + file_name.string() + ".faa";
    }


    std::string rsem(std::string input_path, std::string bam_path, bool paired_end, int threads) {
        // return path
        entapInit::print_msg("Running RSEM...");
        boostFS::path out_dir(ENTAP_EXECUTE::EXECUTE_OUT_PATH+"rsem/");
        boostFS::remove_all(out_dir.c_str());
        boostFS::create_directories(out_dir);
        boostFS::path file_name(input_path);
        file_name = file_name.stem();
        if (file_name.has_stem()) file_name = file_name.stem(); // for .fasta.fnn
        boostFS::path bam_ext(bam_path);bam_ext = bam_ext.extension();
        std::string bam = bam_ext.string();
        std::transform(bam.begin(), bam.end(), bam.begin(), ::tolower);

        if (bam_path.empty()) {
            entapInit::print_msg("No BAM/SAM file provided, exiting RSEM run");
            return input_path;
        }
        if (!entapInit::file_exists(bam_path)) {
            throw ExceptionHandler("Invalid file path for BAM/SAM file, exiting...",
                                   ENTAP_ERR::E_INIT_TAX_READ);
        }
        std::string rsem_arg;
        std::string out_path;
        if (bam.compare(".sam")==0) {
            entapInit::print_msg("File is detected to be sam file, running validation "
                                         "and conversion to bam");
            rsem_arg = ENTAP_EXECUTE::RSEM_EXE_PATH + "rsem-sam-validator " + bam_path;
            out_path = out_dir.string() + file_name.string() + "_rsem_valdate";
            if (entapInit::execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) {
                // only thrown in failure in calling rsem
                throw ExceptionHandler("Error in validating sam file", ENTAP_ERR::E_INIT_TAX_READ);
            }
            // RSEM does not return error code if file is invalid, only seen in .err
            if (!is_file_empty(out_path+".err")) {
                throw ExceptionHandler("Alignment file invalid!", ENTAP_ERR::E_INIT_TAX_READ);
            }
            entapInit::print_msg("Alignment file valid. Converting to BAM");
            std::string bam_out = out_dir.string() + file_name.string();
            rsem_arg = ENTAP_EXECUTE::RSEM_EXE_PATH + "convert-sam-for-rsem " + " -p " + std::to_string(threads)
                       + " "+bam_path + " " + bam_out;
            out_path = out_dir.string() + file_name.string() + "_rsem_convert";
            if (entapInit::execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) {
                // execution error, dif from conversion error
                throw ExceptionHandler("Error in converting sam file", ENTAP_ERR::E_INIT_TAX_READ);
            }
            if (!is_file_empty(out_path+".err")) {
                throw ExceptionHandler("Error in converting sam file", ENTAP_ERR::E_INIT_TAX_READ);
            }
            bam_path = bam_out + ".bam";

        } else if(bam.compare(".bam")==0) {
            entapInit::print_msg("File is detected to be bam file, validating...");
            rsem_arg = ENTAP_EXECUTE::RSEM_EXE_PATH + "rsem-sam-validator " + bam_path;
            out_path = out_dir.string() + file_name.string() + "_rsem_valdate";
            if (entapInit::execute_cmd(rsem_arg.c_str(), out_path.c_str())!=0) {
                throw ExceptionHandler("Error in validating bam file", ENTAP_ERR::E_INIT_TAX_READ);
            }
            if (!is_file_empty(out_path+".err")) {
                throw ExceptionHandler("Alignment file invalid!", ENTAP_ERR::E_INIT_TAX_READ);
            }
            entapInit::print_msg("Alignment file valid. Continuing...");
        } else {
            throw ExceptionHandler("Unknown extension found in the alignment file",
                ENTAP_ERR::E_INIT_TAX_READ);
        }
        // Now have valid BAM file to run rsem
        entapInit::print_msg("Preparing reference");
        std::string ref_out_path = out_dir.string() + file_name.string() + "_ref";
        boostFS::create_directories(ref_out_path);
        rsem_arg = ENTAP_EXECUTE::RSEM_EXE_PATH + "rsem-prepare-reference " + input_path +
                " " + ref_out_path+"/"+file_name.string();
        std::string std_out = out_dir.string() + file_name.string() + "_rsem_reference";
        entapInit::print_msg("Executing following command\n" + rsem_arg);
        entapInit::execute_cmd(rsem_arg.c_str(), std_out);
        entapInit::print_msg("Reference successfully created");

        entapInit::print_msg("Running expression analysis...");
        std::string exp_out_path = out_dir.string() + file_name.string();
        rsem_arg = ENTAP_EXECUTE::RSEM_EXE_PATH + "rsem-calculate-expression " + "--paired-end " +
            "--bam " + "-p " + std::to_string(threads) + " " + bam_path +" "+ ref_out_path +
            "/"+file_name.string()+ " " +exp_out_path;
        if (paired_end) rsem_arg += " --paired-end";
        std_out = out_dir.string() + file_name.string() + "_rsem_exp";
        entapInit::print_msg("Executing following command\n" + rsem_arg);
        if (entapInit::execute_cmd(rsem_arg.c_str(), std_out)!=0) {
            throw ExceptionHandler("Error in running expression analysis",ENTAP_ERR::E_INIT_TAX_READ);
        }
        return out_path + ".genes.results";

    }

    std::string filter_transcriptome(std::string &genemark_path, std::string &rsem_path,
        float fpkm, std::string input_path) {
        entapInit::print_msg("Beginning to filter transcriptome...");
        bool genemark, rsem;

        boostFS::path file_name(input_path);file_name.filename();
        if (genemark_path.empty()) {
            entapInit::print_msg("Looking for genemark file");
            std::string temp_path = ENTAP_EXECUTE::GENEMARK_EXE_PATH+"genemark"+"/"+
                    file_name.string()+ ".faa";
            std::cout<<temp_path<<std::endl;
            if (entapInit::file_exists(temp_path)){
                entapInit::print_msg("File found at: " + temp_path);
                genemark_path = temp_path;
                genemark = true;
            } else {
                entapInit::print_msg("File was not found.");
                genemark = false;
            }
        }else genemark = true;
        if (rsem_path.empty()) {
            entapInit::print_msg("Looking for rsem file");
            std::string temp_path = ENTAP_EXECUTE::EXECUTE_OUT_PATH+"rsem"+"/"+
                                    file_name.stem().string()+ ".genes.results";
            if (entapInit::file_exists(temp_path)){
                entapInit::print_msg("File found at: " + temp_path);
                rsem_path = temp_path;
                rsem = true;
            } else {
                entapInit::print_msg("File was not found.");
                rsem= false;
            }
        }else rsem = true;

        if (!rsem && !genemark) {
            throw ExceptionHandler("Neither genemark, nor rsem files were found", ENTAP_ERR::E_INIT_TAX_READ);
        }
        std::string out_path = ENTAP_EXECUTE::EXECUTE_OUT_PATH + "/" +
                               file_name.stem().string()+"_filtered"+file_name.extension().string();
        std::string out_removed = ENTAP_EXECUTE::EXECUTE_OUT_PATH + "/" +
                                  file_name.stem().string()+"_removed"+file_name.extension().string();
        if (genemark && !rsem) {
            entapInit::print_msg("No rsem file found, so genemark results will continue as main trancriptome: " +
                genemark_path);
            boostFS::copy_file(genemark_path,out_path);
            return genemark_path;
        }

        std::string process_file;
        !genemark ? process_file = input_path: process_file=genemark_path;
        entapInit::print_msg("Filtering file located at: " + process_file);

        io::CSVReader<ENTAP_EXECUTE::RSEM_COL_NUM,io::trim_chars<' ' >,
                io::no_quote_escape<'\t'>> in(rsem_path);
        std::unordered_map<std::string,float> expression_map;
        in.next_line();
        std::string geneid, transid; float length,e_leng, e_count, tpm, fpkm_val;
        while(in.read_row(geneid, transid, length,e_leng,e_count,tpm,fpkm_val)) {
            expression_map.emplace(geneid,fpkm_val);
        }
        remove(out_path.c_str()); remove(out_removed.c_str());
        std::ifstream in_file(process_file);
        std::ofstream out_file(out_path, std::ios::out | std::ios::app);
        std::ofstream removed_file(out_removed, std::ios::out | std::ios::app);
        boost::smatch match;
        bool filtered = false;
        double removed_count = 0;
        for (std::string line; getline(in_file,line);) {
            boost::regex exp(">(\\S+)",boost::regex::icase);
            if (boost::regex_search(line,match,exp)) {
                std::string id = std::string(match[1].first, match[1].second);
                if (expression_map.find(id) != expression_map.end()) {
                    float fp = expression_map.at(id);
                    if (fp > fpkm) {
                        filtered = true;
                        out_file << line;
                    } else  {
                        filtered = false;
                        removed_count++;
                        removed_file << line<<std::endl;
                    }
                } else {
                    // default if not found is NOT remove it
                    filtered = true;
                    out_file << line<<std::endl;
                }
            } else {
                // anything not a seq header
                if (filtered) {
                    out_file << line <<std::endl;
                } else removed_file << line<<std::endl;
            }
        }
        in_file.close(); out_file.close(); removed_file.close();
        entapInit::print_msg("File successfully filtered. Outputs at: " + out_path + " and: " +
            out_removed);
        return out_path;
    }

    void diamond_run(std::list<std::string> database_paths, std::string input_path, int &threads) {
        // not always known (depending on starting state)
        if (!entapInit::file_exists(input_path)) {
            throw ExceptionHandler("Transcriptome file not found",ENTAP_ERR::E_INIT_INDX_DATA_NOT_FOUND);
        }
        boostFS::path transc_name(input_path); transc_name=transc_name.stem();
        if (transc_name.has_stem()) transc_name = transc_name.stem(); //.fasta.faa
        // database verification already ran, don't need to verify each path
        try {
            for (std::string path : database_paths) {
                // assume all paths should be .dmnd
                boostFS::path file_name(path); file_name=file_name.stem();
                entapInit::print_msg("Searching against database located at: " +
                                     path + "...");
                std::string out_path = ENTAP_CONFIG::DIAMOND_RUN_OUT_PATH + transc_name.string() + "_" +
                    file_name.string() + ".out";
                std::string std_out = ENTAP_CONFIG::DIAMOND_RUN_OUT_PATH + transc_name.string() + "_" +
                                      file_name.string() + "_std_";
                diamond_blast(input_path, out_path, std_out, threads);
                entapInit::print_msg("Success! Results written to " + out_path);
            }
        } catch (ExceptionHandler &e) {
            throw ExceptionHandler(e.what(), e.getErr_code());
        }
//        diamond_parse(databases, database_index);
    }

    // input: 3 database string array of selected databases
    void diamond_parse(std::list<std::string> databases, std::string contams) {
        std::unordered_map<std::string, std::string> taxonomic_database;
        std::map<std::string, QuerySequence> query_map;
        try {
            taxonomic_database = read_tax_map();
        } catch (ExceptionHandler &e) {
            throw ExceptionHandler(e.what(),e.getErr_code());
        }
        entapInit::print_msg("Beginning to filter individual databases...");
        for (std::string data : databases) {
            entapInit::print_msg("Database located at "+ data + " being filtered");
            io::CSVReader<ENTAP_EXECUTE::diamond_col_num,io::trim_chars<' ' >,io::no_quote_escape<'\t'>> in(data);
            // todo have columns from input file, in_read_header for versatility
//            in.read_header(io::ignore_extra_column,"qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
//            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle");
            std::string qseqid,sseqid, stitle; float evalue;
            double pident, bitscore; int length, mismatch, gapopen, qstart, qend, sstart, send;
            while(in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                              qstart, qend, sstart ,send, evalue, bitscore, stitle)) {
                QuerySequence new_query = QuerySequence(data, qseqid,sseqid, stitle, evalue);
                boost::regex exp("\\[([^]]+)\\]");      // TODO determined by database format
                boost::smatch match;
                if (boost::regex_search(stitle,match,exp)) {
                    std::string species_lower = std::string(match[1].first, match[1].second);
                    std::string species;strcpy(species,species_lower);
                    std::transform(species_lower.begin(),species_lower.end(), species_lower.begin(),::tolower);
                    std::cout<<"Species found: " +species_lower<<std::endl;
                    if (taxonomic_database.find(species_lower) != taxonomic_database.end()) {
                        new_query.setSpecies(species);
                        new_query.setContaminant(is_contaminant(species_lower, taxonomic_database));
                    } else {
                        new_query.setSpecies("NOT_FOUND");
                        new_query.setContaminant(false);
                    }
                } else{std::cout<<"no match"<<std::endl;}
                // can implement buckets if memory is not issue
                if (query_map.find(qseqid) != query_map.end()) {
                    QuerySequence temp = query_map.at(qseqid);
                    // todo filter database files separately?
                    if (new_query > temp) {
                        query_map.at(qseqid) = new_query;
                    }
                } else {
                    query_map.emplace(qseqid,new_query);
                    std::cout<<"Not found"<<std::endl;
                }
            }
//            print_map(query_map);
        }
    }


    void diamond_blast(std::string input_file, std::string output_file, std::string std_out,int &threads) {
        std::string diamond_run = ENTAP_CONFIG::DIAMOND_PATH_EXE + " blastx " " -d " + input_file +
        " -q " + ENTAP_CONFIG::INPUT_FILE_PATH + " -o " + output_file + " -p " + std::to_string(threads) +" -f " +
                "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle";
        if (entapInit::execute_cmd(diamond_run, std_out) != 0) {
            throw ExceptionHandler("Error in DIAMOND run with database located at: " +
                input_file, ENTAP_ERR::E_INIT_TAX_INDEX);
        }
    }

    std::unordered_map<std::string, std::string> read_tax_map() {
        entapInit::print_msg("Reading taxonomic database into memory...");
        std::unordered_map<std::string, std::string> restored_map;
        try {
            {
                std::ifstream ifs(ENTAP_CONFIG::TAX_BIN_PATH);
                boost::archive::binary_iarchive ia(ifs);
                ia >> restored_map;
            }
        } catch (std::exception &exception){
            throw ExceptionHandler(exception.what(), ENTAP_ERR::E_INIT_TAX_READ);
        }
        entapInit::print_msg("Success!");
        return restored_map;
    }

    bool is_contaminant(std::string species, std::unordered_map<std::string, std::string> &database) {
        // TODO check lineage
        // species and tax database both lowercase



        return false;
    }

    void print_map(std::unordered_map<std::string, QuerySequence> &map) {
        for(std::unordered_map<std::string,QuerySequence>::iterator it = map.begin(); it != map.end(); ++it) {
            QuerySequence q = map.at(it->first);
            std::cout << q;
        }
    }

    void verify_state(std::string &input) {
        std::unordered_map<int, ExecuteStates> enum_map;

        for (int i = 0; input.length(); i++) {
            if (input[i] == 'x') {
            }
            int c = int(input[i]);

        }
    }

    bool is_file_empty(std::string path) {
        std::ifstream ifstream(path);
        return ifstream.peek() == std::ifstream::traits_type::eof();
    }
}