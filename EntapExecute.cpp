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
#include "csv.h"
#include "QuerySequence.h"
#include <boost/regex.hpp>

namespace entapExecute {

    enum ExecuteStates {
        FRAME_SELECTION     = 0x01,
        EXPRESSION          = 0x02,
        DIAMOND_RUN         = 0x04,
        DIAMOND_PARSE       = 0x08,
        EXECUTE_EXIT        = 0x16

    };
    ExecuteStates state = FRAME_SELECTION;

    void execute_main(std::unordered_map<std::string, std::string> user_input) {
        entapInit::print_msg("enTAP Executing...");

        while (state != EXECUTE_EXIT) {
            try {
                genemarkST(user_input.at("i"));
                diamond_run(user_input["U"], user_input["N"], user_input["d"]);

                verify_state(user_input.at("s"));
                state = EXECUTE_EXIT;

            } catch (ExceptionHandler &e) {
                throw ExceptionHandler(e.what(), e.getErr_code());
            }
        }
    }

    void genemarkST(std::string file_path) {

    }
    void rsem(std::string bam) {

    }

    void diamond_run(std::string uniprot, std::string ncbi, std::string database) {
        // Indexed databases
        if (uniprot.empty()) {
            entapInit::print_msg("No Uniprot database selected");
        } else {

        }
        std::string databases[3];
        int database_index = 0;

        std::string uniprot_path = ENTAP_CONFIG::UNIPROT_INDEX_PATH + uniprot + ".dmnd";
        std::string uniprot_out_path = ENTAP_CONFIG::DIAMOND_RUN_OUT_PATH + uniprot + ".out";
        std::string uniprot_std_out_path = ENTAP_CONFIG::DIAMOND_RUN_OUT_PATH + uniprot + "std";

        std::string ncbi_path = ENTAP_CONFIG::NCBI_INDEX_PATH + ncbi + ".dmnd";

        if (!entapInit::file_exists(uniprot_path)) {
            throw ExceptionHandler("Uniprot indexed file not found", ENTAP_ERR::E_INIT_INDX_DATA_NOT_FOUND);
        } else {
            try {
                databases[database_index] = uniprot_out_path;
                database_index++;
                entapInit::print_msg("Searching against Uniprot database located at: " +
                    uniprot_path + "...");
                diamond_blast(uniprot_path, uniprot_out_path, uniprot_std_out_path);
                entapInit::print_msg("Success! Results written to " + uniprot_out_path);
            } catch (ExceptionHandler &e) {
                throw ExceptionHandler(e.what(), e.getErr_code());
            }
        }
        diamond_parse(databases, database_index);

    }

    // input: 3 database string array of selected databases
    void diamond_parse(std::list databases) {
        std::unordered_map<std::string, std::string> taxonomic_database;
        std::unordered_map<std::string, QuerySequence> query_map;
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
                boost::regex exp("\\[([^]]+)\\]");
                boost::smatch match;
                if (boost::regex_search(stitle,match,exp)) {
                    std::string species = std::string(match[1].first, match[1].second);
                    std::transform(species.begin(),species.end(), species.begin(),::tolower);
                    std::cout<<"Species found: " +species<<std::endl;
                    if (taxonomic_database.find(species) != taxonomic_database.end()) {
                        new_query.setSpecies(species);
                        new_query.setContaminant(is_contaminant(species));
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
            print_map(query_map);
        }
    }


    void diamond_blast(std::string input_file, std::string output_file, std::string std_out) {
        std::string diamond_run = ENTAP_CONFIG::DIAMOND_PATH_EXE + " blastx " " -d " + input_file +
        " -q " + ENTAP_CONFIG::INPUT_FILE_PATH + " -o " + output_file +" -f " +
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

    bool is_contaminant(std::string cont) {
        return false;
    }

    void print_map(std::unordered_map<std::string, QuerySequence> map) {
        for(std::unordered_map<std::string,QuerySequence>::iterator it = map.begin(); it != map.end(); ++it) {
            QuerySequence q = map.at(it->first);
            std::cout << q;
        }
    }

    void verify_state(std::string input) {
        std::unordered_map<int, ExecuteStates> enum_map;

    }

}