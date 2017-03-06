//
// Created by harta on 3/4/17.
//

#include <boost/serialization/unordered_map.hpp>
#include <iostream>
#include "EntapExecute.h"
#include "ExceptionHandler.h"
#include "EntapConsts.h"
#include "EntapInit.h"

namespace entapExecute {

    enum ExecuteStates {

    };

    void execute_main(std::unordered_map<std::string, std::string> user_input) {
        std::cout << "EXECUTE" << std::endl;
        try {
            diamond_run(user_input["U"], user_input["N"], user_input["d"]);
            diamond_parse();

        } catch (ExceptionHandler &e) {
            throw ExceptionHandler(e.what(), e.getErr_code());
        }
    }

    void genemarkST() {

    }
    void rsem() {

    }

    void diamond_run(std::string uniprot, std::string ncbi, std::string database) {
        std::string uniprot_path = ENTAP_CONST::UNIPROT_INDEX_PATH + uniprot + ".dmnd";
        std::string uniprot_out_path = ENTAP_CONST::DIAMOND_RUN_OUT_PATH + uniprot + ".out";

        std::string ncbi_path = ENTAP_CONST::NCBI_INDEX_PATH + ncbi + ".dmnd";

        if (!entapInit::file_exists(uniprot_path)) {
            throw ExceptionHandler("Uniprot indexed file not found", ENTAP_ERR::E_INIT_INDX_DATA_NOT_FOUND);
        } else {
            std::string diamond_uniprot_run = ENTAP_CONST::DIAMOND_PATH_EXE + " blastx " " -d " + uniprot_path +
                " -q " + ENTAP_CONST::INPUT_FILE_PATH + " -o " + uniprot_out_path;
            entapInit::execute_cmd(diamond_uniprot_run);

        }
    }

    void diamond_parse() {

    }

}