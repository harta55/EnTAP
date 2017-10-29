/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017, Alexander Hart, Dr. Jill Wegrzyn
 *
 * This file is part of EnTAP.
 *
 * EnTAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EnTAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ENTAPCONFIG_H
#define ENTAPCONFIG_H

#include <boost/serialization/unordered_map.hpp>
#include <boost/program_options/variables_map.hpp>
#include "common.h"

struct  GoEntry {
    std::string go_id;
    std::string level;
    std::string category;
    std::string term;
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const uint32 v) {
        ar&go_id;
        ar&level;
        ar&category;
        ar&term;
    }
};

struct TaxEntry {
    std::string tax_id;
    std::string lineage;
    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive & ar, const uint32 v) {
        ar&tax_id;
        ar&lineage;
    }
    bool is_empty() {
        return this->tax_id.empty() && this->lineage.empty();
    }
};

namespace entapConfig {
    //-----------------------FTP PATHS---------------------------//
    const std::string GO_DATABASE_FTP =
            "http://archive.geneontology.org/latest-full/go_monthly-termdb-tables.tar.gz";
    const std::string UNIPROT_FTP_SWISS = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";
    const std::string UNIPROT_FTP_TREMBL = "ftp://ftp.uniprot.org/pub/databases/"
            "uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz";
    const std::string GO_TERM_FILE = "term.txt";
    const std::string GO_GRAPH_FILE = "graph_path.txt";
    const std::string GO_DATA_NAME = "go_monthly-termdb-tables.tar.gz";
    const std::string GO_DIR = "go_monthly-termdb-tables/";
    const std::string ENTAP_CONFIG_DIR = "/entap_config";
    const std::string TAX_DATABASE_PATH = "/databases/ncbi_tax.entp";

    //******************Prototype Functions******************

    void init_entap(boost::program_options::variables_map, std::string);
    void init_taxonomic(std::string&);
    void init_uniprot(std::vector<std::string>&, std::string);
    void init_ncbi(std::vector<std::string>&, std::string);
    void init_diamond_index(std::string,std::string,int);
    std::string download_file(std::string, std::string&,std::string&);
    std::string download_file(const std::string &,std::string&);
    void decompress_file(std::string,std::string,short);
    int update_database(std::string);
    void init_go_db(std::string&,std::string);
    void init_eggnog(std::string);
    void handle_state(void);
}

#endif //ENTAP_INITHANDLER_H
