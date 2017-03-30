//
// Created by harta on 3/4/17.
//

#include <boost/serialization/unordered_map.hpp>
#include <iostream>
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapConsts.h"
#include "EntapInit.h"
#include "csv.h"

namespace entapExecute {

    enum ExecuteStates {

    };

    void execute_main(std::unordered_map<std::string, std::string> user_input) {
        entapInit::print_msg("enTAP Executing...");
        try {
            diamond_run(user_input["U"], user_input["N"], user_input["d"]);

        } catch (ExceptionHandler &e) {
            throw ExceptionHandler(e.what(), e.getErr_code());
        }
    }

    void genemarkST() {

    }
    void rsem() {

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
    void diamond_parse(std::string databases[3], int s) {
        entapInit::print_msg("Beginning to filter individual databases...");
        for (int i = 0; i < s; i++) {
            entapInit::print_msg("Database located at "+ databases[i] + " being filtered");
            io::CSVReader<ENTAP_EXECUTE::diamond_col_num> in(databases[i]);
            std::string line[13];
            std::string str[50];
            int q = 0;
            while(in.read_row(line)) {
                str[q] = line[0];
                q++;
            }
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

}