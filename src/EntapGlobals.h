/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
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
#include <boost/serialization/access.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "config.h"
#include "common.h"

//******************** Defines/Macros **************************

#ifdef USE_BOOST
#define PATHS(x,y)      (boostFS::path(x) / boostFS::path(y)).string()
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


namespace ENTAP_EXECUTE {
    //------------------------Ontology-------------------------//
    extern const std::string GO_BIOLOGICAL_FLAG ;
    extern const std::string GO_CELLULAR_FLAG;
    extern const std::string GO_MOLECULAR_FLAG;
    const uint16 ONTOLOGY_MIN         = 0;
#ifdef EGGNOG_MAPPER
    const uint16 EGGNOG_INT_FLAG      = 0;
#endif
    const uint16 INTERPRO_INT_FLAG    = 1;
    const uint16 EGGNOG_DMND_INT_FLAG = 0;  // Set to 0 for EggNOG / mapper
    const uint16 ONTOLOGY_MAX         = 1;

    const uint16 FRAME_FLAG_GENEMARK = 0;
    const uint16 EXP_FLAG_RSEM       = 0;
    const uint16 SIM_SEARCH_FLAG_DIAMOND = 0;

    const uint16 SOFTWARE_MAX = 2;


    //------------------------Headers-------------------------//

    // Sim Search Header Information
    extern const std::string HEADER_QUERY;
    extern const std::string HEADER_SUBJECT;
    extern const std::string HEADER_PERCENT;
    extern const std::string HEADER_ALIGN_LEN;
    extern const std::string HEADER_MISMATCH;
    extern const std::string HEADER_GAP_OPEN;
    extern const std::string HEADER_QUERY_S;
    extern const std::string HEADER_QUERY_E;
    extern const std::string HEADER_SUBJ_S;
    extern const std::string HEADER_SUBJ_E;
    extern const std::string HEADER_E_VAL;
    extern const std::string HEADER_COVERAGE;
    extern const std::string HEADER_TITLE;
    extern const std::string HEADER_SPECIES;
    extern const std::string HEADER_DATABASE;
    extern const std::string HEADER_FRAME;
    extern const std::string HEADER_CONTAM;
    extern const std::string HEADER_INFORM;

    // UniProt Mapping Header Information
    extern const std::string HEADER_UNI_DATA_XREF;
    extern const std::string HEADER_UNI_COMMENTS;

    // EggNOG Header Information
    extern const std::string HEADER_SEED_ORTH;
    extern const std::string HEADER_SEED_EVAL;
    extern const std::string HEADER_SEED_SCORE;
    extern const std::string HEADER_PRED_GENE;
    extern const std::string HEADER_TAX_SCOPE;
    extern const std::string HEADER_EGG_OGS;
    extern const std::string HEADER_EGG_KEGG;
    extern const std::string HEADER_EGG_GO_BIO ;
    extern const std::string HEADER_EGG_GO_CELL;
    extern const std::string HEADER_EGG_GO_MOLE;
    extern const std::string HEADER_EGG_DESC;
    extern const std::string HEADER_EGG_LEVEL;
    extern const std::string HEADER_EGG_PROTEIN;

    // Interpro Header Information
    extern const std::string HEADER_INTER_GO_BIO;
    extern const std::string HEADER_INTER_GO_CELL;
    extern const std::string HEADER_INTER_GO_MOLE;
    extern const std::string HEADER_INTER_PATHWAY;
    extern const std::string HEADER_INTER_INTERPRO;
    extern const std::string HEADER_INTER_DATA_TYPE;
    extern const std::string HEADER_INTER_DATA_TERM;
    extern const std::string HEADER_INTER_EVAL;
}

namespace ENTAP_STATS {
    extern const std::string SOFTWARE_BREAK;

    void ES_format_stat_stream(std::stringstream &stream, std::string title);


}

#endif //ENTAPGLOBALS_H
