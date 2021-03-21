/******************************************************************
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
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
 *******************************************************************/

#ifndef ENTAP_BUSCODATABASE_H
#define ENTAP_BUSCODATABASE_H

#include "../common.h"

class FileSystem;

/**
 * ======================================================================
 * @class BuscoDatabase
 *
 * Description          - This module provides an interface to the BUSCO
 *                        database which provides for an assessment of the
 *                        transcriptome/annotation utilizing OrthoDB
 *                      - The user may input various datasets that are
 *                        analyzed through BUSCO
 *
 * Dependencies         - Python (python3 is preferred according to documentation)
 *                      - BLAST
 *
 * Citation             - BUSCO applications from quality assessments to gene
 *                        prediction and phylogenomics. Robert M. Waterhouse,
 *                        Mathieu Seppey, Felipe A. Simão, Mose Manni, Panagiotis
 *                        Ioannidis, Guennadi Klioutchnikov, Evgenia V. Kriventseva,
 *                        and Evgeny M. Zdobnov Mol Biol Evol, published online
 *                        Dec 6, 2017 | doi: 10.1093/molbev/msx319
 *
 * ======================================================================
 */
class BuscoDatabase {

public:

    typedef enum {

        ERR_DATA_OK,
        ERR_DATA_UNKNOWN_DATABASE,
        ERR_DATA_DOWNLOAD,
        ERR_DATA_DECOMPRESS


    } BUSCO_DB_ERR;

    BuscoDatabase(FileSystem* fileSystem);
    ~BuscoDatabase();

    BUSCO_DB_ERR download_database(std::string& user_input, std::string& output_path);
    bool valid_database(std::string &database, std::string &database_url);

    // Error handling
    std::string print_error_log();

private:
    std::string get_database_shortname(std::string &database);

    // Error handling
    void set_err_msg(std::string msg, BUSCO_DB_ERR code);

    // URL mappings taken from https://busco.ezlab.org/frame_wget.html
    const std::unordered_map<std::string, std::string> DATASET_TO_URL = {

            // Bacteria URLs
            {"bacteria"                 , "http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz"},
            {"proteobacteria"           , "http://busco.ezlab.org/v2/datasets/proteobacteria_odb9.tar.gz"},
            {"rhizobiales"              , "http://busco.ezlab.org/v2/datasets/rhizobiales_odb9.tar.gz"},
            {"betaproteobacteria"       , "http://busco.ezlab.org/v2/datasets/betaproteobacteria_odb9.tar.gz"},
            {"gammaproteobacteria", "http://busco.ezlab.org/v2/datasets/gammaproteobacteria_odb9.tar.gz"},
            {"enterobacteriales", "http://busco.ezlab.org/v2/datasets/enterobacteriales_odb9.tar.gz"},
            {"deltaepsilonsub", "http://busco.ezlab.org/v2/datasets/deltaepsilonsub_odb9.tar.gz"},
            {"actinobacteria" , "http://busco.ezlab.org/v2/datasets/actinobacteria_odb9.tar.gz"},
            {"cyanobacteria", "http://busco.ezlab.org/v2/datasets/cyanobacteria_odb9.tar.gz"},
            {"firmicutes",    "http://busco.ezlab.org/v2/datasets/firmicutes_odb9.tar.gz" },
            {"clostridia",      "http://busco.ezlab.org/v2/datasets/clostridia_odb9.tar.gz"},
            {"lactobacillales", "http://busco.ezlab.org/v2/datasets/lactobacillales_odb9.tar.gz"},
            {"bacillales", "http://busco.ezlab.org/v2/datasets/bacillales_odb9.tar.gz"},
            {"bacteroidetes", "http://busco.ezlab.org/v2/datasets/bacteroidetes_odb9.tar.gz"},
            {"spirochaetes", "http://busco.ezlab.org/v2/datasets/spirochaetes_odb9.tar.gz"},
            {"tenericutes", "http://busco.ezlab.org/v2/datasets/tenericutes_odb9.tar.gz"},

            // Eukaryota URLs
            {"eukaryota", "http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz"},
            {"fungi", "http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz"},
            {"microsporidia", "http://busco.ezlab.org/v2/datasets/microsporidia_odb9.tar.gz"},
            {"dikarya", "http://busco.ezlab.org/v2/datasets/dikarya_odb9.tar.gz"},
            {"ascomycota", "http://busco.ezlab.org/v2/datasets/ascomycota_odb9.tar.gz"},
            {"pezizomycotina", "http://busco.ezlab.org/v2/datasets/pezizomycotina_odb9.tar.gz"},
            {"eurotiomycetes", "http://busco.ezlab.org/v2/datasets/eurotiomycetes_odb9.tar.gz"},
            {"sordariomyceta", "http://busco.ezlab.org/v2/datasets/sordariomyceta_odb9.tar.gz"},
            {"saccharomyceta", "http://busco.ezlab.org/v2/datasets/saccharomyceta_odb9.tar.gz"},
            {"saccharomycetales", "http://busco.ezlab.org/v2/datasets/saccharomycetales_odb9.tar.gz"},
            {"basidiomycota", "http://busco.ezlab.org/v2/datasets/basidiomycota_odb9.tar.gz"},
            {"metazoa", "http://busco.ezlab.org/v2/datasets/metazoa_odb9.tar.gz"},
            {"nematoda", "http://busco.ezlab.org/v2/datasets/nematoda_odb9.tar.gz"},
            {"arthropoda", "http://busco.ezlab.org/v2/datasets/arthropoda_odb9.tar.gz"},
            {"insecta", "http://busco.ezlab.org/v2/datasets/insecta_odb9.tar.gz"},
            {"endopterygota", "http://busco.ezlab.org/v2/datasets/endopterygota_odb9.tar.gz"},
            {"hymenoptera", "http://busco.ezlab.org/v2/datasets/hymenoptera_odb9.tar.gz"},
            {"diptera", "http://busco.ezlab.org/v2/datasets/diptera_odb9.tar.gz"},
            {"vertebrata", "http://busco.ezlab.org/v2/datasets/vertebrata_odb9.tar.gz"},
            {"actinopterygii", "http://busco.ezlab.org/v2/datasets/actinopterygii_odb9.tar.gz"},
            {"tetrapoda", "http://busco.ezlab.org/v2/datasets/tetrapoda_odb9.tar.gz"},
            {"aves", "http://busco.ezlab.org/v2/datasets/aves_odb9.tar.gz"},
            {"mammalia", "http://busco.ezlab.org/v2/datasets/mammalia_odb9.tar.gz"},
            {"euarchontoglires", "http://busco.ezlab.org/v2/datasets/euarchontoglires_odb9.tar.gz"},
            {"laurasiatheria", "http://busco.ezlab.org/v2/datasets/laurasiatheria_odb9.tar.gz"},
            {"embryophyta", "http://busco.ezlab.org/v2/datasets/embryophyta_odb9.tar.gz"},
            {"protists", "http://busco.ezlab.org/v2/datasets/protists_ensembl.tar.gz"},
            {"alveolata_stramenophiles", "http://busco.ezlab.org/v2/datasets/alveolata_stramenophiles_ensembl.tar.gz"}
    };


    FileSystem *mpFileSystem;
    std::string mTempDirectory;     // Absolute path to temp directory to download files
    std::string mErrMsg;
    BUSCO_DB_ERR mErrCode;
};


#endif //ENTAP_BUSCODATABASE_H
