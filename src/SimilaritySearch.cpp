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

//*********************** Includes *****************************
#include <boost/range/iterator_range_core.hpp>
#include <csv.h>
#include <boost/regex.hpp>
#include "SimilaritySearch.h"
//**************************************************************

typedef std::map<std::string,std::map<std::string,uint32>> graph_sum_t;


/**
 * ======================================================================
 * Function SimilaritySearch(std::vector<std::string> &databases, std::string input,
                           int threads, std::string out, boost::program_options::variables_map &user_flags,
                           GraphingManager *graphingManager)
 *
 * Description          - SimilaritySearch object constructor
 *                      - Responsible for initiating SimSearch member variables,
 *                        parsing user input for relevant information used within
 *                        this module
 *
 * Notes                - None
 *
 * @param databases     - List of user selected databases (parsed in  execute namespace)
 * @param input         - Path to user transcriptome (final version post filtering)
 * @param threads       - Thread count
 * @param out           - EnTAP out directory
 * @param user_flags    - Boost parsed user input
 * @param graphingManager- Pointer to graphing manager
 *
 * @return              - SimilaritySearch instance
 * ======================================================================
 */
SimilaritySearch::SimilaritySearch(databases_t &databases,
                                   std::string input, EntapDataPtrs& entap_data) {
    FS_dprint("Spawn Object - SimilaritySearch");
    std::string uninform_path;

    _pQUERY_DATA    = entap_data._pQueryData;
    _pUserInput     = entap_data._pUserInput;
    _pFileSystem    = entap_data._pFileSystem;
    _pEntapDatabase = entap_data._pEntapDatbase;
    _pGraphingManager = entap_data._pGraphingManager;
    _database_paths = databases;
    _input_path     = input;
    _diamond_exe    = DIAMOND_EXE;      // Set to extern set previously

    // Species already checked for validity in Init
    _input_species    = _pUserInput->get_target_species_str();
    _qcoverage        = _pUserInput->get_user_input<fp32>(UInput::INPUT_FLAG_QCOVERAGE);
    _tcoverage        = _pUserInput->get_user_input<fp32>(UInput::INPUT_FLAG_TCOVERAGE);
    _overwrite        = _pUserInput->has_input(UInput::INPUT_FLAG_OVERWRITE);
    _e_val            = _pUserInput->get_user_input<fp64>(UInput::INPUT_FLAG_E_VAL);
    _threads          = _pUserInput->get_supported_threads();
    _uninformative_vect= _pUserInput->get_uninformative_vect();
    _outpath          = _pFileSystem->get_root_path();
    _contaminants     = _pUserInput->get_contaminants();
    _software_flag    = ENTAP_EXECUTE::SIM_SEARCH_FLAG_DIAMOND; // Default DIAMOND software

    // Set sim search paths/directories
    _sim_search_dir  = PATHS(_outpath, SIM_SEARCH_DIR);
    _processed_path  = PATHS(_sim_search_dir, PROCESSED_DIR);
    _results_path    = PATHS(_sim_search_dir, RESULTS_DIR);

    if (_overwrite) {
        _pFileSystem->delete_dir(_sim_search_dir);
    }

    _pFileSystem->create_dir(_sim_search_dir);
    _pFileSystem->delete_dir(_processed_path);
    _pFileSystem->create_dir(_processed_path);
    _pFileSystem->delete_dir(_results_path);
    _pFileSystem->create_dir(_results_path);

    _transcript_shortname = get_transcriptome_shortname();
}


/**
 * ======================================================================
 * Function SimilaritySearch()
 *
 * Description          - Empty constructor
 *
 * Notes                - None
 *
 * @return              - SimilaritySearch instance
 * ======================================================================
 */
SimilaritySearch::SimilaritySearch() {}


/**
 * ======================================================================
 * Function std::vector<std::string> SimilaritySearch::execute(std::string updated_input,
 *                                                             bool blast)
 *
 * Description          - Responsible for executing the selected similarity
 *                        searching software
 *                      - Returns output files of sim search
 *                      - Throws instance of ExceptionHandler on failed execution
 *
 * Notes                - None
 *
 * @param updated_input - User transcriptome, if changed
 * @param blast         - Blastx/blastp (true for blastp)
 *
 * @return              - Vector of output files
 * ======================================================================
 */
std::vector<std::string> SimilaritySearch::execute(std::string updated_input,bool blast) {
    this->_input_path = updated_input;
    this->_blastp = blast;
    _blastp ? _blast_type = BLASTP : _blast_type = BLASTX;
    try {
        switch (_software_flag) {
            case 0:
                return diamond();
            default:
                return diamond();
        }
    } catch (ExceptionHandler &e) {throw e;}
}


/**
 * ======================================================================
 * Function std::pair<std::string,std::string> SimilaritySearch::parse_files(std::string new_input,
                                   std::map<std::string, QuerySequence>& MAP)
 *
 * Description          - Responsible for best hit selection of sequences hit
 *                        against the diamond databases
 *                      - Throws fatal ExceptionHandler instance
 *
 * Notes                - Entered as SIM_SEARCH_PARSE state from Execute namespace
 *
 * @param new_input     - User transcriptome, if changed
 * @param MAP           - Master data structure of transcriptome information
 *
 * @return              - Pair of output files (hits and no hits)
 * ======================================================================
 */
void SimilaritySearch::parse_files(std::string new_input) {
    _input_path = new_input;
    try {
        switch (_software_flag) {
            case ENTAP_EXECUTE::SIM_SEARCH_FLAG_DIAMOND:
                diamond_parse(_contaminants);
                break;
            default:
                diamond_parse(_contaminants);
                break;
        }
    } catch (ExceptionHandler &e) {throw e;}
}


/**
 * ======================================================================
 * Function std::vector<std::string> SimilaritySearch::diamond()
 *
 * Description          - Responsible for executing simliarity search through
 *                        pstreams library
 *                      - Returns vector of output files from sim search
 *                      - Checks whether DIAMOND has been ran previously
 *
 * Notes                - None
 *
 *
 * @return              - Output files from similarity searching
 * ======================================================================
 */
std::vector<std::string> SimilaritySearch::diamond() {
    FS_dprint("Beginning to execute DIAMOND...");

    std::vector<std::string>    out_paths;
    std::string                 filename;
    std::string                 out_path;
    std::string                 std_out;
    std::string                 database_name;  // shortened name

    if (!_pFileSystem->file_exists(_input_path)) {
        throw ExceptionHandler("Transcriptome file not found",ERR_ENTAP_RUN_SIM_SEARCH_RUN);
    }

    // database verification already ran, don't need to verify each path
    try {
        // assume all paths should be .dmnd
        for (std::string data_path : _database_paths) {
            FS_dprint("Searching against database located at: " + data_path + "...");
            database_name = get_database_shortname(data_path);
            filename = _blast_type + "_" + _transcript_shortname + "_" + database_name + FileSystem::EXT_OUT;
            out_path = PATHS(_sim_search_dir,filename) ;
            std_out  = PATHS(_sim_search_dir,filename) + "_std";
            _file_to_database[out_path] = database_name;
            if (_pFileSystem->file_exists(out_path)) {
                FS_dprint("File found at " + out_path + " skipping execution against this database");
                out_paths.push_back(out_path);
                continue;
            }
            diamond_blast(_input_path, out_path, std_out,data_path, _threads, _blast_type);
            FS_dprint("Success! Results written to " + out_path);
            out_paths.push_back(out_path);
        }
    } catch (const ExceptionHandler &e) {throw e;}
    _sim_search_paths = out_paths;
    return out_paths;
}


/**
 * ======================================================================
 * Function void SimilaritySearch::diamond_blast(std::string input_file, std::string output_file, std::string std_out,
                   std::string &database,int &threads, std::string &blast)
 *
 * Description          - Responsible for execution of DIAMOND through pstreams
 *                        library
 *
 * Notes                - None
 *
 * @param input_file    - Path to input transcriptome
 * @param output_file   - Path to output file from sim search
 * @param std_out       - Std out/err path
 * @param database      - Selected database to hit against
 * @param threads       - Thread number
 * @param blast         - Blast type (blastx/blastp)
 *
 * @return              - None
 * ======================================================================
 */
void SimilaritySearch::diamond_blast(std::string input_file, std::string output_file, std::string std_out,
                   std::string &database,int &threads, std::string &blast) {

    std::string        diamond_run;
    std::stringstream  stream_err;
    std::stringstream  stream_out;
    int32              err_code;

    diamond_run =
            _diamond_exe + " "
            + blast +
            " -d " + database    +
            " --query-cover "    + std::to_string(_qcoverage) +
            " --subject-cover "  + std::to_string(_tcoverage) +
            " --evalue "         + std::to_string(_e_val) +
            " --more-sensitive"  +
            " --top 3"           +
            " -q " + input_file  +
            " -o " + output_file +
            " -p " + std::to_string(threads) +
            " -f " + "6 qseqid sseqid pident length mismatch gapopen "
                     "qstart qend sstart send evalue bitscore qcovhsp stitle";

    // Execute command
    err_code = TC_execute_cmd(diamond_run, stream_err, stream_out, std_out);
    FS_dprint("\nDIAMOND STD OUT:\n" + stream_out.str());
    if (err_code != 0) {
        // Delete output file if run failed
        _pFileSystem->delete_file(output_file);
        throw ExceptionHandler("Error with database located at: " + database + "\n" +
                                       stream_err.str(), ERR_ENTAP_RUN_SIM_SEARCH_RUN);
    }
}


/**
 * ======================================================================
 * Function std::vector<std::string> SimilaritySearch::verify_diamond_files(
 *                                      std::string &outpath,
 *                                      std::string name)
 *
 * Description          - Checks whether DIAMOND has already been ran with
 *                        the same parameters to skip re-running
 *
 * Notes                - None
 *
 * @param outpath       - Path to out directory
 * @param name          - Name of input transcriptome
 *
 * @return              - None
 * ======================================================================
 */
std::vector<std::string> SimilaritySearch::verify_diamond_files() {
    FS_dprint("Override unselected, checking for diamond files of selected databases...");
    std::vector<std::string> out_list;
    std::string              temp_out;
    std::string              file_name_full;
    std::string              database_name;     // shortened version

    for (std::string &data_path : _database_paths) {
        // assume all paths should be .dmnd
        database_name = get_database_shortname(data_path);
        file_name_full = _blast_type + "_" + _transcript_shortname + "_" + database_name + FileSystem::EXT_OUT;
        temp_out = PATHS(_sim_search_dir, file_name_full);
        if (!_pFileSystem->file_exists(temp_out)){
            FS_dprint("File at: " + temp_out + " not found, running diamond");
            out_list.clear();
            return out_list;
        }
        out_list.push_back(temp_out);
    }
    FS_dprint("All diamond files found, skipping this stage of EnTAP");
    _sim_search_paths = out_list;
    return out_list;
}

// input: 3 database string array of selected databases
void SimilaritySearch::diamond_parse(std::vector<std::string>& contams) {
    FS_dprint("Beginning to filter individual diamond_files...");

    std::list<std::map<std::string,QuerySequence>>  database_maps;
    std::pair<bool,std::string>                     contam_info;
    std::string                                     species;
    TaxEntry                                        taxEntry;
    QuerySequence::SimSearchResults                 simSearchResults;
    bool                                            is_uniprot;
    uint32                                          uniprot_attempts=0;

    // ------------------ Read from DIAMOND output ---------------------- //
    std::string qseqid;
    std::string sseqid, stitle, database_name,pident, bitscore,
            length, mismatch, gapopen, qstart, qend, sstart, send;
    fp64  evalue;
    fp64  coverage;

    // ----------------------------------------------------------------- //

    // Get the taxonomic info (lineage) of the target species
    taxEntry = _pEntapDatabase->get_tax_entry(_input_species);
    _input_lineage = taxEntry.lineage;

    for (std::string &data : _sim_search_paths) {
        // Loop through each database
        FS_dprint("Diamond file located at " + data + " being filtered");
        std::stringstream out_stream;
        is_uniprot = false;
        uniprot_attempts = 0;

        // Confirm we have legit path / not empty
        if (!_pFileSystem->file_exists(data) || _pFileSystem->file_empty(data)) {
            // Should never fall into here
            throw ExceptionHandler("File not found or empty: " + data, ERR_ENTAP_RUN_SIM_SEARCH_FILTER);
        }

        // Set up individual database directories for stats
        database_name = _file_to_database[data];

        // Begin using CSVReader lib to parse data
        io::CSVReader<DMND_COL_NUMBER, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(data);
        while (in.read_row(qseqid, sseqid, pident, length, mismatch, gapopen,
                           qstart, qend, sstart, send, evalue, bitscore, coverage,stitle)) {
            simSearchResults = {};

            // Get pointer to sequence in overall map
            QuerySequence *query = _pQUERY_DATA->get_sequence(qseqid);
            if (query == nullptr) {
                throw ExceptionHandler("Unable to find sequence in transcriptome: " + qseqid + " from file: " + data,
                                       ERR_ENTAP_RUN_SIM_SEARCH_FILTER);
            }

            // get species from database alignment (using boost regex for now)
            species = get_species(stitle);
            // get taxonomic information with species
            taxEntry = _pEntapDatabase->get_tax_entry(species);
            // get contaminant information
            contam_info = is_contaminant(taxEntry.lineage,contams);

            // Check if this is a UniProt match and pull back info if so
            if (is_uniprot) {
                // Get uniprot info
                is_uniprot_entry(sseqid, simSearchResults.uniprot_info);
            } else {
                if (uniprot_attempts <= UNIPROT_ATTEMPTS) {
                    // First UniProt match assumes the rest are UniProt as well
                    is_uniprot = is_uniprot_entry(sseqid, simSearchResults.uniprot_info);
                    if (!is_uniprot) uniprot_attempts++;
                }
            }


            // Compile sim search data
            simSearchResults.database_path = data;
            simSearchResults.qseqid = qseqid;
            simSearchResults.sseqid = sseqid;
            simSearchResults.pident = pident;
            simSearchResults.length = length;
            simSearchResults.mismatch = mismatch;
            simSearchResults.gapopen = gapopen;
            simSearchResults.qstart = qstart;
            simSearchResults.qend = qend;
            simSearchResults.sstart = sstart;
            simSearchResults.send = send;
            simSearchResults.stitle = stitle;
            simSearchResults.bit_score = bitscore;
            simSearchResults.lineage = taxEntry.lineage;
            simSearchResults.species = species;
            simSearchResults.e_val_raw = evalue;
            simSearchResults.e_val = float_to_sci(evalue,2);
            simSearchResults.coverage_raw = coverage;
            simSearchResults.coverage = float_to_string(coverage);
            simSearchResults.contaminant = contam_info.first;
            simSearchResults.contam_type = contam_info.second;
            simSearchResults.contaminant ? simSearchResults.yes_no_contam = "Yes" :
                    simSearchResults.yes_no_contam  = "No";
            simSearchResults.is_informative = is_informative(stitle);
            simSearchResults.is_informative ? simSearchResults.yes_no_inform = "Yes" :
                    simSearchResults.yes_no_inform  = "No";

            query->add_alignment(
                    SIMILARITY_SEARCH,
                    _software_flag,
                    simSearchResults,
                    data,
                    _input_lineage);
        }

        FS_dprint("File parsed, calculating statistics and writing output...");
        calculate_best_stats(false,data);
        FS_dprint("Success!");
    }
    FS_dprint("Calculating overall Similarity Searching statistics...");
    calculate_best_stats(true);
    FS_dprint("Success!");
}

void SimilaritySearch::calculate_best_stats (bool is_final, std::string database_path) {

    GraphingData                graphingStruct;
    std::string                 species;
    std::string                 database_shortname;
    std::string                 figure_base;
    std::string                 frame;
    std::string                 base_path;
    std::string                 temp_file_path;
    std::string                 contam;
    std::stringstream           ss;
    uint64                      count_no_hit=0;
    uint64                      count_contam=0;
    uint64                      count_filtered=0;
    uint64                      count_informative=0;
    uint64                      count_uninformative=0;
    uint64                      count_unselected=0;
    uint64                      count_TOTAL_alignments=0;
    uint32                      ct;
    fp64                        percent;
    fp64                        contam_percent;
    Compair<std::string>        contam_counter;
    Compair<std::string>        species_counter;
    Compair<std::string>        contam_species_counter;
    graph_sum_t                 graphing_sum_map;

    // Set up output directories (processed directory cleared earlier so these will be empty)
    if (is_final) {
        // Overall results across databases
        base_path = _results_path;
        database_shortname = "";
    } else {
        // Individual database results
        database_shortname = _file_to_database[database_path];
        base_path   = PATHS(_processed_path, database_shortname);
    }
    figure_base = PATHS(base_path, FIGURE_DIR);
    _pFileSystem->create_dir(base_path);
    _pFileSystem->create_dir(figure_base);

    // Open contam best hit tsv file
    std::string out_best_contams_tsv = PATHS(base_path ,SIM_SEARCH_DATABASE_CONTAM_TSV);
    std::ofstream file_best_contam_tsv(out_best_contams_tsv,std::ios::out | std::ios::app);

    // Open contam fasta protein file
    std::string out_best_contams_fa_prot = PATHS(base_path ,SIM_SEARCH_DATABASE_CONTAM_FA_PROT);
    std::ofstream file_best_contam_fa_prot(out_best_contams_fa_prot,std::ios::out | std::ios::app);

    // Open contam fasta nucleotide file
    temp_file_path           = PATHS(base_path ,SIM_SEARCH_DATABASE_CONTAM_FA_NUCL);
    std::ofstream file_best_contam_fa_nucl(temp_file_path,std::ios::out | std::ios::app);

    // Open best hits file (tsv)
    std::string out_best_hits_tsv         = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_TSV);
    std::ofstream file_best_hits_tsv(out_best_hits_tsv,std::ios::out | std::ios::app);

    // Open best hits file (fasta nucleotide)
    std::string out_best_hits_fa_nucl = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_FA_NUCL);
    std::ofstream file_best_hits_fa_nucl(out_best_hits_fa_nucl,std::ios::out | std::ios::app);

    // Open best hits file (fasta protein)
    std::string out_best_hits_fa_prot    = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_FA_PROT);
    std::ofstream file_best_hits_fa_prot(out_best_hits_fa_prot,std::ios::out | std::ios::app);

    // Open best hits file with no contaminants (tsv)
    temp_file_path           = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_TSV_NO_CONTAM);
    std::ofstream file_best_hits_tsv_no_contam(temp_file_path,std::ios::out | std::ios::app);

    // Open best hits file with no contaminants (fasta nucleotide)
    temp_file_path           = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_FA_NUCL_NO_CONTAM);
    std::ofstream file_best_hits_fa_nucl_no_contam(temp_file_path,std::ios::out | std::ios::app);

    // Open best hits file with no contaminants (fasta protein)
    temp_file_path           = PATHS(base_path, SIM_SEARCH_DATABASE_BEST_FA_PROT_NO_CONTAM);
    std::ofstream file_best_hits_fa_prot_no_contam(temp_file_path,std::ios::out | std::ios::app);

    // Open unselected hits, so every hit that was not the best hit (tsv)
    std::string out_unselected_tsv  = PATHS(base_path, SIM_SEARCH_DATABASE_UNSELECTED);
    std::ofstream file_unselected_hits(out_unselected_tsv, std::ios::out | std::ios::app);

    // Open no hits file (fasta nucleotide)
    std::string out_no_hits_fa_nucl = PATHS(base_path, SIM_SEARCH_DATABASE_NO_HITS_NUCL);
    std::ofstream file_no_hits_nucl(out_no_hits_fa_nucl, std::ios::out | std::ios::app);

    // Open no hits file (fasta protein)
    std::string out_no_hits_fa_prot  = PATHS(base_path, SIM_SEARCH_DATABASE_NO_HITS_PROT);
    std::ofstream file_no_hits_prot(out_no_hits_fa_prot, std::ios::out | std::ios::app);

    // ------------------- Setup graphing files ------------------------- //

    std::string graph_species_txt_path           = PATHS(figure_base, GRAPH_SPECIES_BAR_TXT);
    std::string graph_species_png_path           = PATHS(figure_base, GRAPH_SPECIES_BAR_PNG);
    std::string graph_contam_txt_path            = PATHS(figure_base, GRAPH_CONTAM_BAR_TXT);
    std::string graph_contam_png_path            = PATHS(figure_base, GRAPH_CONTAM_BAR_PNG);
    std::string graph_sum_txt_path               = PATHS(figure_base, GRAPH_DATABASE_SUM_TXT);
    std::string graph_sum_png_path               = PATHS(figure_base, GRAPH_DATABASE_SUM_PNG);

    std::ofstream graph_species_file(graph_species_txt_path, std::ios::out | std::ios::app);
    std::ofstream graph_contam_file(graph_contam_txt_path, std::ios::out | std::ios::app);
    std::ofstream graph_sum_file(graph_sum_txt_path, std::ios::out | std::ios::app);

    // ------------------------------------------------------------------ //

    // Print headers to relevant tsv files
    print_header(file_best_contam_tsv);
    print_header(file_best_hits_tsv);
    print_header(file_best_hits_tsv_no_contam);
    print_header(file_unselected_hits);


    try {
        graph_species_file << "Species\tCount"     << std::endl;
        graph_contam_file  << "Contaminant Species\tCount" << std::endl;
        graph_sum_file     << "Category\tCount"    << std::endl;

        // Cycle through all sequences
        for (auto &pair : *_pQUERY_DATA->get_sequences_ptr()) {
            // Check if original sequences have hit a database
            if (!pair.second->hit_database(SIMILARITY_SEARCH, ENTAP_EXECUTE::SIM_SEARCH_FLAG_DIAMOND, database_path)) {
                // Did NOT hit a database during sim search
                // Do NOT log if it was never blasted
                if ((pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN) && _blastp) ||
                        (!pair.second->QUERY_FLAG_GET(QuerySequence::QUERY_IS_PROTEIN) && !_blastp)) {
                    // Protein/nucleotide did not hit database
                    count_no_hit++;
                    file_no_hits_nucl << pair.second->get_sequence_n() << std::endl;
                    file_no_hits_prot << pair.second->get_sequence_p() << std::endl;
                    // Graphing
                    frame = pair.second->getFrame();
                    if (graphing_sum_map[frame].find(NO_HIT_FLAG) != graphing_sum_map[frame].end()) {
                        graphing_sum_map[frame][NO_HIT_FLAG]++;
                    } else graphing_sum_map[frame][NO_HIT_FLAG] = 1;
                } else {
                    pair.second->QUERY_FLAG_SET(QuerySequence::QUERY_BLASTED);
                }
            } else {
                // HIT a database during sim search

                QuerySequence::SimSearchResults *sim_search_data;
                QuerySequence::SimSearchAlignment *best_hit;
                // Process unselected hits for non-final analysis and set best hit pointer
                if (is_final) {
                    best_hit =
                            pair.second->get_best_hit_alignment<QuerySequence::SimSearchAlignment>(
                                    SIMILARITY_SEARCH,ENTAP_EXECUTE::SIM_SEARCH_FLAG_DIAMOND,"");
                    sim_search_data = best_hit->get_results();
                } else {
                    best_hit = pair.second->get_best_hit_alignment<QuerySequence::SimSearchAlignment>(
                            SIMILARITY_SEARCH,ENTAP_EXECUTE::SIM_SEARCH_FLAG_DIAMOND,database_path);
                    QuerySequence::align_database_hits_t *alignment_data =
                            pair.second->get_database_hits(database_path,SIMILARITY_SEARCH, ENTAP_EXECUTE::SIM_SEARCH_FLAG_DIAMOND);
                    sim_search_data = best_hit->get_results();
                    for (auto &hit : alignment_data->second) {
                        count_TOTAL_alignments++;
                        if (hit != best_hit) {  // If this hit is not the best hit
                            file_unselected_hits << hit->print_tsv(DEFAULT_HEADERS) << std::endl;
                            count_unselected++;
                        } else {
                            ;   // Do notthing
                        }
                    }
                }
                count_filtered++;   // increment best hit
                // Write to best hits files
                file_best_hits_fa_nucl << pair.second->get_sequence_n()<<std::endl;
                file_best_hits_fa_prot << pair.second->get_sequence_p()<<std::endl;
                file_best_hits_tsv << best_hit->print_tsv(DEFAULT_HEADERS) << std::endl;

                frame = pair.second->getFrame();     // Used for graphing
                species = sim_search_data->species;

                // Determine contaminant information and print to files
                if (sim_search_data->contaminant) {
                    // Species is considered a contaminant
                    count_contam++;
                    file_best_contam_fa_nucl << pair.second->get_sequence_n()<<std::endl;
                    file_best_contam_fa_prot << pair.second->get_sequence_p()<<std::endl;
                    file_best_contam_tsv << best_hit->print_tsv(DEFAULT_HEADERS) << std::endl;
                    contam = sim_search_data->contam_type;
                    contam_counter.add_value(contam);
                    contam_species_counter.add_value(species);
                } else {
                    // Species is NOT a contaminant, print to files
                    file_best_hits_fa_nucl_no_contam << pair.second->get_sequence_n()<<std::endl;
                    file_best_hits_fa_prot_no_contam << pair.second->get_sequence_p()<<std::endl;
                    file_best_hits_tsv_no_contam << best_hit->print_tsv(DEFAULT_HEADERS) << std::endl;
                }

                // Count species type
                species_counter.add_value(species);

                // Check if this is an informative alignment and respond accordingly
                if (sim_search_data->is_informative) {
                    count_informative++;
                    // Graphing
                    if (graphing_sum_map[frame].find(INFORMATIVE_FLAG) != graphing_sum_map[frame].end()) {
                        graphing_sum_map[frame][INFORMATIVE_FLAG]++;
                    } else graphing_sum_map[frame][INFORMATIVE_FLAG] = 1;

                } else {
                    count_uninformative++;
                    if (graphing_sum_map[frame].find(UNINFORMATIVE_FLAG) != graphing_sum_map[frame].end()) {
                        graphing_sum_map[frame][UNINFORMATIVE_FLAG]++;
                    } else graphing_sum_map[frame][UNINFORMATIVE_FLAG] = 1;
                }
            }
        }
    } catch (const std::exception &e){throw ExceptionHandler(e.what(), ERR_ENTAP_RUN_SIM_SEARCH_FILTER);}

    try {
        _pFileSystem->close_file(file_best_hits_tsv);
        _pFileSystem->close_file(file_best_hits_tsv_no_contam);
        _pFileSystem->close_file(file_best_hits_fa_nucl);
        _pFileSystem->close_file(file_best_hits_fa_prot);
        _pFileSystem->close_file(file_best_hits_fa_nucl_no_contam);
        _pFileSystem->close_file(file_best_hits_fa_prot_no_contam);
        _pFileSystem->close_file(file_best_contam_tsv);
        _pFileSystem->close_file(file_best_contam_fa_prot);
        _pFileSystem->close_file(file_best_contam_fa_nucl);
        _pFileSystem->close_file(file_no_hits_nucl);
        _pFileSystem->close_file(file_no_hits_prot);
        _pFileSystem->close_file(file_unselected_hits);
    } catch (const ExceptionHandler &e) {throw e;}

    // ------------ Calculate statistics and print to output ------------ //
    ss<<std::fixed<<std::setprecision(2);

    // Different headers if final analysis or database specific analysis
    if (is_final) {
        ss << ENTAP_STATS::SOFTWARE_BREAK
           << "Compiled Similarity Search - Diamond - Best Overall" << "\n"
           << ENTAP_STATS::SOFTWARE_BREAK;
    } else {
        ss << ENTAP_STATS::SOFTWARE_BREAK
           << "Similarity Search - Diamond - " << database_shortname << "\n"
           << ENTAP_STATS::SOFTWARE_BREAK <<
           "Search results:\n"            << database_path <<
           "\n\tTotal alignments: "               << count_TOTAL_alignments   <<
           "\n\tTotal unselected results: "       << count_unselected      <<
           "\n\t\tWritten to: "                   << out_unselected_tsv;
    }

    // If overall alignments are 0, then throw error
    if (is_final && count_filtered == 0) {
        throw ExceptionHandler("No alignments found during Similarity Searching!",
                               ERR_ENTAP_RUN_SIM_SEARCH_FILTER);
    }

    // If no total or filealignments for this database, return and warn user
    if (!is_final && (count_TOTAL_alignments == 0 || count_filtered == 0)) {
        ss << "WARNING: No alignments for this database";
        std::string out_msg = ss.str() + "\n";
        _pFileSystem->print_stats(out_msg);
        return;
    }

    // Sort counters
    contam_species_counter.sort(true);
    species_counter.sort(true);

    contam_percent = ((fp64) count_contam / count_filtered) * 100;

    ss <<
       "\n\tTotal unique transcripts with an alignment: " << count_filtered <<
       "\n\t\tReference transcriptome sequences with an alignment (FASTA):\n\t\t\t" << out_best_hits_fa_prot <<
       "\n\t\tSearch results (TSV):\n\t\t\t" << out_best_hits_tsv <<
       "\n\tTotal unique transcripts without an alignment: " << count_no_hit <<
       "\n\t\tReference transcriptome sequences without an alignment (FASTA):\n\t\t\t" << out_no_hits_fa_prot;
    // Have frame information
    if (graphing_sum_map.size() > 1) {
        for (auto &pair : graphing_sum_map) {
            // Frame -> Map of uninform/inform/no hits
            ss << "\n\t\t" << pair.first << "(" << pair.second[NO_HIT_FLAG] << ")";
            graph_sum_file << pair.first << "\t" << NO_HIT_FLAG << "\t" << pair.second[NO_HIT_FLAG] << "\n";
        }
    }
    ss <<
       "\n\tTotal unique informative alignments: " << count_informative;
    if (graphing_sum_map.size() > 1) {
        for (auto &pair : graphing_sum_map) {
            // Frame -> Map of uninform/inform/no hits
            ss << "\n\t\t" << pair.first << "(" << pair.second[INFORMATIVE_FLAG] << ")";
            graph_sum_file << pair.first << "\t" << INFORMATIVE_FLAG << "\t" << pair.second[INFORMATIVE_FLAG]
                           << "\n";
        }
    }
    ss <<
       "\n\tTotal unique uninformative alignments: " << count_uninformative;
    if (graphing_sum_map.size() > 1) {
        for (auto &pair : graphing_sum_map) {
            // Frame -> Map of uninform/inform/no hits
            ss << "\n\t\t" << pair.first << "(" << pair.second[UNINFORMATIVE_FLAG] << ")";
            graph_sum_file << pair.first << "\t" << UNINFORMATIVE_FLAG << "\t" << pair.second[UNINFORMATIVE_FLAG]
                           << "\n";
        }
    }

    ss <<
       "\n\tTotal unique contaminants: " << count_contam <<
       "(" << contam_percent << "%): " <<
       "\n\t\tTranscriptome reference sequences labeled as a contaminant (FASTA):\n\t\t\t"
       << out_best_contams_fa_prot <<
       "\n\t\tTranscriptome reference sequences labeled as a contaminant (TSV):\n\t\t\t" << out_best_contams_tsv;


    // ********** Contaminant Calculations ************** //
    if (count_contam > 0) {
        ss << "\n\t\tFlagged contaminants (all % based on total contaminants):";
        for (auto &pair : contam_counter._data) {
            percent = ((fp64) pair.second / count_contam) * 100;
            ss
                    << "\n\t\t\t" << pair.first << ": " << pair.second << "(" << percent << "%)";
        }
        ss << "\n\t\tTop " << COUNT_TOP_SPECIES << " contaminants by species:";
        ct = 1;
        for (auto &pair : contam_species_counter._sorted) {
            if (ct > COUNT_TOP_SPECIES) break;
            percent = ((fp64) pair.second / count_contam) * 100;
            ss
                    << "\n\t\t\t" << ct << ")" << pair.first << ": "
                    << pair.second << "(" << percent << "%)";
            graph_contam_file << pair.first << '\t' << std::to_string(pair.second) << std::endl;
            ct++;
        }
    }

    ss << "\n\tTop " << COUNT_TOP_SPECIES << " alignments by species:";
    ct = 1;
    for (auto &pair : species_counter._sorted) {
        if (ct > COUNT_TOP_SPECIES) break;
        percent = ((fp64) pair.second / count_filtered) * 100;
        ss
                << "\n\t\t\t" << ct << ")" << pair.first << ": "
                << pair.second << "(" << percent << "%)";
        graph_species_file << pair.first << '\t' << std::to_string(pair.second) << std::endl;
        ct++;
    }
    std::string out_msg = ss.str() + "\n";
    _pFileSystem->print_stats(out_msg);


    // ------------------------------------------------------------------ //
    // ********* Graphing Handle ********** //
    graphingStruct.software_flag = GRAPH_SOFTWARE_FLAG;
    _pFileSystem->close_file(graph_contam_file);
    _pFileSystem->close_file(graph_species_file);
    _pFileSystem->close_file(graph_sum_file);
    if (count_contam > 0) {
        graphingStruct.fig_out_path   = graph_contam_png_path;
        graphingStruct.graph_title    = database_shortname + GRAPH_CONTAM_TITLE;
        graphingStruct.text_file_path = graph_contam_txt_path;
        graphingStruct.graph_type     = GRAPH_BAR_FLAG;
        _pGraphingManager->graph(graphingStruct);
    }
    graphingStruct.fig_out_path   = graph_species_png_path;
    graphingStruct.graph_title    = database_shortname + GRAPH_SPECIES_TITLE;
    graphingStruct.text_file_path = graph_species_txt_path;
    graphingStruct.graph_type     = GRAPH_BAR_FLAG;
    _pGraphingManager->graph(graphingStruct);

    graphingStruct.fig_out_path   = graph_sum_png_path;
    graphingStruct.graph_title    = database_shortname + GRAPH_DATABASE_SUM_TITLE;
    graphingStruct.text_file_path = graph_sum_txt_path;
    graphingStruct.graph_type     = GRAPH_SUM_FLAG;
    _pGraphingManager->graph(graphingStruct);

    // check if final - different graph
    // ************************************ //
}

std::pair<bool,std::string> SimilaritySearch::is_contaminant(std::string lineage,
                    std::vector<std::string> &contams) {
    // species and tax database both lowercase
    if (contams.empty()) return std::pair<bool,std::string>(false,"");
    LOWERCASE(lineage);
    for (auto const &contaminant:contams) {
        if (lineage.find(contaminant) != std::string::npos){
            return std::pair<bool,std::string>(true,contaminant);
        }
    }
    return std::pair<bool,std::string>(false,"");
}

std::string SimilaritySearch::get_species(std::string &title) {
    // TODO use regex(database specific)

    std::string species;
    boost::smatch match;

    boost::regex ncbi_exp(_NCBI_REGEX);
    boost::regex uniprot_exp(_UNIPROT_REGEX);

    species = "";
    if (boost::regex_search(title,match,uniprot_exp)) {
        species = std::string(match[1].first, match[1].second);
    } else {
        if (boost::regex_search(title, match, ncbi_exp))
            species = std::string(match[1].first, match[1].second);
    }
    // Double bracket fix
    if (species[0] == '[') species = species.substr(1);
    if (species[species.length()-1] == ']') species = species.substr(0,species.length()-1);
    return species;
}

bool SimilaritySearch::is_informative(std::string title) {
    LOWERCASE(title);
    for (std::string &item : _uninformative_vect) { // Already lowercase
        if (title.find(item) != std::string::npos) return false;
    }
    return true;
}

void SimilaritySearch::print_header(std::ofstream &file_stream) {
    for (const std::string *header : DEFAULT_HEADERS) {
        file_stream << *header << "\t";
    }
    file_stream << std::endl;
}

bool SimilaritySearch::is_executable() {
    std::string test_command;

    test_command = DIAMOND_EXE + " --version";
    return TC_execute_cmd(test_command) == 0;
}

// Returns database "shortname" from database full path
std::string SimilaritySearch::get_database_shortname(std::string &path) {
    return boostFS::path(path).filename().stem().string();
}

std::string SimilaritySearch::get_transcriptome_shortname() {
    boostFS::path transc_name(_input_path);
    transc_name = transc_name.stem();
    if (transc_name.has_stem()) transc_name = transc_name.stem(); //.fasta.faa
    return transc_name.string();
}

bool SimilaritySearch::is_uniprot_entry(std::string &sseqid, UniprotEntry &entry) {
    std::string accession;

    // sseqid - sp|Q9FJZ9|PER72_ARATH
    if (_pEntapDatabase == nullptr || sseqid.empty()) return false;

    accession = sseqid.substr(sseqid.rfind('|',sseqid.length())+1);     // Q9FJZ9
    entry = _pEntapDatabase->get_uniprot_entry(accession);
    return !entry.is_empty();
}

