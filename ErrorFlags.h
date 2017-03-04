//
// Created by harta55 on 2/1/17.
//

#ifndef ENTAP_ERRORFLAGS_H
#define ENTAP_ERRORFLAGS_H

namespace ENTAPERR {
    const int E_INPUT_PARSE = 10;
    const int E_SUCCESS = 11;
    const int E_INIT_TAX_DOWN = 20;
    const int E_INIT_TAX_INDEX = 21;
    const int E_INIT_TAX_SERIAL = 22;
    const int E_INIT_TAX_READ = 55;


    const std::string taxonomic_script = "../download_tax.pl";
    const std::string taxonomic_database = "databases/ncbi_tax.entp";
    const std::string UNIPROT_DEFAULT = "swiss";
    const std::string NCBI_DEFAULT = "uniref";

    const std::string UNIPROT_FTP_SWISS = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";
    con

}



#endif //ENTAP_ERRORFLAGS_H
