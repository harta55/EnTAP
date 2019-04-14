/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2019, Alexander Hart, Dr. Jill Wegrzyn
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


#ifndef ENTAPGLOBALS_H
#define ENTAPGLOBALS_H

//*********************** Includes *****************************

#include "config.h"
#include "common.h"
#ifdef USE_BOOST
#include <boost/serialization/access.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#endif
//******************** Defines/Macros **************************

#ifdef USE_BOOST
#define PATHS(x,y)      (boostFS::path(x) / boostFS::path(y)).string()
#else
#define PATHS(x,y)      ((x) + "/" + (y))
#endif
#define FASTA_FLAG      ">"

//**************************************************************


class QuerySequence;
class EntapDatabase;
class FileSystem;
class UserInput;
class GraphingManager;
class QueryData;

#ifdef USE_BOOST
namespace boostFS = boost::filesystem;
namespace boostPO = boost::program_options;
namespace boostAR = boost::archive;
#endif
//**************************************************************



//***************** Global Prototype Functions *****************
std::string generate_command(std::unordered_map<std::string,std::string>&,
                             std::string);
std::string float_to_string(fp64);
std::string float_to_sci(fp64, int);
vect_str_t  split_string(std::string, char);
std::string get_cur_time();
//**************************************************************


//**************** Global Structures/Typedefs ******************
typedef std::unordered_map<std::string, QuerySequence*> QUERY_MAP_T;
typedef std::map<std::string,std::vector<std::string>> go_format_t;
typedef std::vector<std::string> databases_t;   // Standard database container


template <typename T>
struct Compair {

    typedef std::pair<T,uint32> count_pair;
    std::unordered_map<T, uint32> _data;
    typename std::unordered_map<T, uint32>::iterator _it;
    std::vector<count_pair> _sorted;
    uint32 _ct_unique;
    uint32 _ct_total;

    struct sort_descending{
        bool operator ()(count_pair const& one, count_pair const& two) const {
            return one.second > two.second;
        }
    };

    struct sort_ascending {
        bool operator ()(count_pair const& one, count_pair const& two) const {
            return one.second < two.second;
        }
    };

    Compair() {
        _data = {};
        _sorted = {};
        _ct_unique = 0;
        _ct_total = 0;
    }

    bool add_value(T val) {
        _ct_total++;
        _it = _data.find(val);
        if (_it != _data.end()) {
            _it->second++;
            return true;
        } else {
            _data.emplace(std::make_pair(val, 1));
            _ct_unique++;
            return false;
        }
    }

    void sort(bool descending) {
        _sorted = std::vector<count_pair>(_data.begin(), _data.end());
        if (descending) {
            std::sort(_sorted.begin(), _sorted.end(), sort_descending());
        } else {
            std::sort(_sorted.begin(), _sorted.end(), sort_ascending());
        }
    }

    bool empty() {
        return _data.empty();
    }
};

template <class T>
std::string container_to_string(std::set<T> &in_cont, const char *delim) {
    if (in_cont.empty()) return "";
    std::ostringstream stream;
    std::copy(in_cont.begin(), in_cont.end(), std::ostream_iterator<std::string>(stream, delim));
    std::string result = stream.str();
    result.pop_back();
    return result;
}

template <template<class,class,class...> class C, typename T, typename U, typename... Args>
U get_map_default(const C<T,U,Args...>& m, T const& key, const U & default_val)
{
    typename C<T,U,Args...>::const_iterator it = m.find(key);
    if (it == m.end()) {
        return default_val;
    }
    return it->second;
}


enum ExecuteStates {
    INIT = 0,
    EXPRESSION_FILTERING,
    FRAME_SELECTION,
    FILTER,
    SIMILARITY_SEARCH,
    GENE_ONTOLOGY,
    EXIT,
    EXECUTION_MAX
};

enum ONTOLOGY_SOFTWARE {
    ONT_EGGNOG_DMND,
    ONT_INTERPRO_SCAN,
#ifdef EGGNOG_MAPPER
    EGGNOG_INT_FLAG,
#endif
    ONT_SOFTWARE_COUNT
};

enum SIMILARITY_SOFTWARE {
    SIM_DIAMOND,
    SIM_SOFTWARE_COUNT
};

enum EXPRESSION_SOFTWARE {
    EXP_RSEM,
    EXP_COUNT
};

enum FRAME_SELECTION_SOFTWARE {
    FRAME_GENEMARK_ST,
    FRAME_SOFTWARE_COUNT
};

enum ENTAP_HEADERS {
    ENTAP_HEADER_UNUSED = 0,                // 0
    ENTAP_HEADER_QUERY,

    /* Frame Selection */
    ENTAP_HEADER_FRAME,

    /* Expression Filtering */
    ENTAP_HEADER_EXP_FPKM,

    /* Similarity Search - General */
    ENTAP_HEADER_SIM_SUBJECT,
    ENTAP_HEADER_SIM_PERCENT,       // Percent Identical to subject
    ENTAP_HEADER_SIM_ALIGN_LEN,
    ENTAP_HEADER_SIM_MISMATCH,
    ENTAP_HEADER_SIM_GAP_OPEN,
    ENTAP_HEADER_SIM_QUERY_S,
    ENTAP_HEADER_SIM_QUERY_E,               // 10
    ENTAP_HEADER_SIM_SUBJ_S,
    ENTAP_HEADER_SIM_SUBJ_E,
    ENTAP_HEADER_SIM_E_VAL,
    ENTAP_HEADER_SIM_COVERAGE,
    ENTAP_HEADER_SIM_TITLE,
    ENTAP_HEADER_SIM_SPECIES,
    ENTAP_HEADER_SIM_DATABASE,
    ENTAP_HEADER_SIM_CONTAM,
    ENTAP_HEADER_SIM_INFORM,

    /* Similarity Search - UniProt*/
    ENTAP_HEADER_SIM_UNI_DATA_XREF,         // 20
    ENTAP_HEADER_SIM_UNI_COMMENTS,
    ENTAP_HEADER_SIM_UNI_KEGG,
    ENTAP_HEADER_SIM_UNI_GO_BIO,
    ENTAP_HEADER_SIM_UNI_GO_CELL,
    ENTAP_HEADER_SIM_UNI_GO_MOLE,

    /* Ontology - EggNOG*/
    ENTAP_HEADER_ONT_EGG_SEED_ORTHO,
    ENTAP_HEADER_ONT_EGG_SEED_EVAL,
    ENTAP_HEADER_ONT_EGG_SEED_SCORE,
    ENTAP_HEADER_ONT_EGG_PRED_GENE,
    ENTAP_HEADER_ONT_EGG_TAX_SCOPE,         // 30
    ENTAP_HEADER_ONT_EGG_OGS,
    ENTAP_HEADER_ONT_EGG_DESC,
    ENTAP_HEADER_ONT_EGG_KEGG,
    ENTAP_HEADER_ONT_EGG_GO_BIO,
    ENTAP_HEADER_ONT_EGG_GO_CELL,
    ENTAP_HEADER_ONT_EGG_GO_MOLE,
    ENTAP_HEADER_ONT_EGG_PROTEIN,

    /* Ontology - InterProScan */
    ENTAP_HEADER_ONT_INTER_GO_BIO,
    ENTAP_HEADER_ONT_INTER_GO_CELL,
    ENTAP_HEADER_ONT_INTER_GO_MOLE,
    ENTAP_HEADER_ONT_INTER_PATHWAYS,
    ENTAP_HEADER_ONT_INTER_INTERPRO,
    ENTAP_HEADER_ONT_INTER_DATA_TYPE,
    ENTAP_HEADER_ONT_INTER_DATA_TERM,
    ENTAP_HEADER_ONT_INTER_EVAL,

    ENTAP_HEADER_COUNT
};

struct EntapHeader {
    const std::string title;
    bool print_header;
};

struct EntapDataPtrs {
    EntapDatabase* _pEntapDatbase;
    FileSystem*    _pFileSystem;
    UserInput*     _pUserInput;
    GraphingManager* _pGraphingManager;
    QueryData*     _pQueryData;

    bool is_null() {
        return _pEntapDatbase == nullptr || _pFileSystem == nullptr ||
        _pUserInput == nullptr || _pGraphingManager == nullptr ||
        _pQueryData == nullptr;
    }

    EntapDataPtrs() {
        _pEntapDatbase = nullptr;
        _pFileSystem   = nullptr;
        _pUserInput    = nullptr;
        _pGraphingManager = nullptr;
        _pQueryData    = nullptr;
    }
};

namespace std {
    template<>
    struct hash<ENTAP_HEADERS> {
        std::size_t operator()(const ENTAP_HEADERS &x) const {
            using type = typename std::underlying_type<ENTAP_HEADERS>::type;
            return std::hash<type>()(static_cast<type>(x));
        }
    };
}



//*********************** Externs *****************************

extern std::string DEBUG_FILE_PATH;
extern std::string LOG_FILE_PATH;
extern std::string RSEM_EXE_DIR;
extern std::string GENEMARK_EXE;
extern std::string DIAMOND_EXE;
extern std::string EGG_SQL_DB_PATH;
extern std::string EGG_DMND_PATH;
extern std::string INTERPRO_EXE;
extern std::string ENTAP_DATABASE_BIN_PATH;
extern std::string ENTAP_DATABASE_SQL_PATH;
extern std::string GRAPHING_EXE;
extern std::string EGG_EMAPPER_EXE;

extern EntapHeader ENTAP_HEADER_INFO[];

// ************************************************************

static const std::string GO_MOLECULAR_FLAG     = "molecular_function";
static const std::string GO_BIOLOGICAL_FLAG    = "biological_process";
static const std::string GO_CELLULAR_FLAG      = "cellular_component";
static const std::string GO_OVERALL_FLAG       = "overall";

#endif //ENTAPGLOBALS_H
