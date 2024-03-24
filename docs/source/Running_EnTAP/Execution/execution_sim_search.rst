.. |sim_dir| replace:: :file:`/similarity_search/DIAMOND`
.. |sim_proc_dir| replace:: :file:`/similarity_search/DIAMOND/processed/database_ref`
.. |sim_overall_dir| replace:: :file:`/similarity_search/DIAMOND/overall_results`
.. |sim_fig_dir| replace:: :file:`/similarity_search/DIAMOND/processed/database_ref/figures`
.. |sim_res_dir| replace:: */overall_results*
.. |ncbi_refseq| replace:: https://www.ncbi.nlm.nih.gov/refseq/
.. |uniprot_swiss| replace:: https://www.uniprot.org/
.. |diamond_git| replace:: https://github.com/bbuchfink/diamond
.. |ncbi_tax| replace:: https://www.ncbi.nlm.nih.gov/taxonomy

Similarity Searching
=========================
Similarity Searching is the process of aligning the sequences at this stage of the pipeline (after filtering and refinement has been performed) to reference databases. These databases often come from NCBI RefSeq (|ncbi_refseq| or UniProt Swiss-Prot (|uniprot_swiss|). The assumption here is that the user has already configured some of these databases for use with DIAMOND. If not, these can be done during the Configuration stage of EnTAP, refer back to that page for more information.

DIAMOND (|diamond_git|) is utilized for Similarity Searching over traditional methods, such as BLAST, due to its improved speeds and accuracy. 

Running Similarity Searching
----------------------------------
Similarity Searching is a mandatory stage of EnTAP. In order to run this, the user must input multiple DIAMOND configured databases (should have the '.dmnd' extension) through the :file:`database` command. There are many customizable parameters for this stage as well that can be seen in the table below. An important piece to consider, is that EnTAP will utilize a unique method of :ref:`best alignment selection<best_hit-label>` that selects the best aligment per database utilizing a variety of parameters. This also factors in contaminant filtering!

Similarity Search Commands
-------------------------------------

.. list-table:: **Similarity Search Flags**
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1    
   
   * - param
     - description
     - location (cmd/R-ini,E-ini)
     - qualifier
     - example
   * - database / d
     - Specify up to 5 DIAMOND indexed (.dmnd) databases to run similarity search against
     - R-ini
     - multi-string
     - /path/to/diamond/database.dmnd
   * - data-type
     - Specify which EnTAP database you'd like to use for execution
         * 0. Binary Database (default) - This will be much quicker and is recommended
         * 1. SQL Database - Slower although will be more easily compatible with every system
     - R-ini
     - integer
     - 0
   * - contam / c
     - Specify contaminants to be used during Simlilarity Search best hit selection. Contaminants can be selected by species or through a specific taxon (insecta) from the NCBI Taxonomy Database. If your taxon is more than one word just replace the spaces with underscores (_). Alignments will be flagged as contaminants and will be lower scoring compared to other alignments.
     - R-ini
     - multi-string
     - insecta
   * - taxon
     - This flag will allow for 'taxonomic favoring' of hits that are closer to your target species or lineage. Any lineage can be used as referenced by the NCBI Taxonomic database, such as genus, phylum, or species. Format **must** replace all spaces with underscores ('_')
     - R-ini
     - string
     - homo_sapiens
   * - e
     - Specify E-value cutoff for Similarity Searching results (in scientific notation format).
     - R-ini
     - scientific
     - 10E-5
   * - tcoverage
     - Specify minimum target coverage for similarity searching
     - R-ini
     - float
     - 50
   * - qcoverage
     - Specify minimum query coverage for similarity searching
     - R-ini
     - float
     - 50
   * - uninformative
     - Comma-deliminated list of terms you would like to be deemed "uninformative". Any alignments during Similarity Searching tagged as uninformative will be scored lower
     - R-ini
     - string
     - conserved, predicted, unnamed, hypothetical, putative, unidentified, uncharacterized, unknown, uncultured, uninformative
   * - diamond-exe
     - Specify the execution method for DIAMOND. This can be a path to the :file:`diamond` file generated during installation, or simply the command if installed globally
     - E-ini
     - string
     - diamond


.. _best_hit-label:

EnTAP Best Alignment Selection
-----------------------------------
EnTAP incorporates a unique method of selecting the best alignment per database, then overall. This utilizes parameters such as E-Value, Coverage, :ref:`Contaminant Filtering<contam-label>`, and :ref:`Taxonomic Favoring<taxon_favor-label>`. 

Best Alignment Selection For Each Database:

    #. Examine E-Value
	
        * If E-Value difference is high, select the smallest E-Value alignment
		
        * If E-Value difference is low, continue
		
    #. Examine Coverage
	
        * If Coverage difference is high, select the larger Coverage alignment
		
        * If Coverage difference is low, continue
		
    #. :ref:`Contaminant Filtering<contam-label>`
	
        * If one alignment is a contaminant, select non-contaminant
		
        * If both/neither alignments are contaminants, continue
		
    #. :ref:`Taxonomic Favoring<taxon_favor-label>` and :ref:`Informativeness<inform-label>`
	
        * Select alignment that is closer in lineage to our target species and more informative

.. _contam-label:

Contaminant Filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Contaminant Filtering leverages the NCBI Taxonomic Database (|ncbi_tax|) to tag alignments that could be considered a contaminant. Contaminants can be introduced during collection or processing of a sample. A contaminant is essentially a species that is not of the target species you are collecting. Some common contaminants are bacteria and fungi that can sometimes be found within collected samples. Oftentimes, researchers would like to remove these sequences from the dataset. 

In order to use Contaminant Filtering, the user must use the :file:`contam` flag to select multiple (comma-separated) species/taxon from the NCBI Taxonomy Database. When inputting these, any spaces in the taxon must be replaced by an underscore ('_'). Any alignments containing the contaminant taxon will be flagged as such and unfavored during best alignment selection. 

As an example with common contaminants within the |run_ini_file_format| file:

.. code-block:: bash
    
    contam=insecta,fungi,bacteria

.. note:: Sometimes the best alignment can be a contaminant! EnTAP will flag this and allow the user to decide whether or not they would like to retain it


.. _taxon_favor-label:

Taxonomic Favoring
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
During best alignment selection of similarity searched results, taxonomic consideration can utilized. If a certain lineage (pinus) is specified, hits closer in taxonomic lineage to this selection will be chosen. Any lineage such as species/kingdom/phylum can be utilized as long as it is contained within the NCBI Taxonomic Database. If it is not located within the database, EnTAP will stop the execution immediately and let you know! 

This feature can be utilized via the |run_ini_file_format| file. An example can be seen below (remember to replace any spaces with an underscore):

.. code-block:: bash

    taxon=pinus_taeda

.. _inform-label:

Informativeness
^^^^^^^^^^^^^^^^^^^^^^^^^
Informativeness is another metric that is used during selection of the best alignment. Oftentimes reference databases may have terms such as 'unknown' in the descriptions of alignments where certain information may not be known about that alignment. EnTAP will attempt to select alignments that are well established rather than these 'uninformative' alignments.

Any term can be used via the :file:`uninformative` flag in the |run_ini_file_format|, so you are not limited! Below are defaults used by EnTAP as comma-separated:

.. code-block:: bash

    uninformative=conserved, predicted, unnamed, hypothetical, putative, unidentified, uncharacterized, unknown, uncultured, uninformative

Interpreting the Results
-------------------------------
The |sim_dir| directory will contain all of the relevant information for the Similarity Searching stage of the pipeline. This folder will contain the files generated from DIAMOND as well as files generated from EnTAP. Files generated from EnTAP for each individual reference database are contained within the |sim_proc_dir|, while the overall analysis compiling the results of each reference database are contained within the |sim_overall_dir|. 

The same files are repeated across databases and across the overall results, so I will only go into detail for each file once below with an input transcriptome labelled 'species' and a reference database labelled 'ref_database'.

.. list-table:: **Similarity Search Results**
   :align: left
   :widths: 10 50 10
   :header-rows: 1    
   
   * - filename
     - description
     - directory
   * - :file:`blastp_species_ref_database.out`
     - Generated from DIAMOND. Contains a lot of information from the DIAMOND search including e-value, coverage, reference database descriptions, and much more. This is a typical output file from a BLAST type of search. The filename is prepended with either blastp or blastx depending on if runP or runN was used. A file like this will be generated for each reference database used.
     - |sim_dir|
   * - :file:`blastp_species_ref_database_std.err/.out`
     - Generated from DIAMOND. These files are will contain any error or general information produced from the DIAMOND run.
     - |sim_dir|
   * - :file:`diamond_annotated.faa/.fnn/.tsv`
     - Generated from EnTAP. Contains all of the best alignments (protein/nucleotide format) that were selected from this database, or overall combining the results from each database used. Since this contains all best alignments, it may contain contaminants or uninformative alignments. Sometimes a contaminant can be the best alignment! Note: Protein or nucleotide information may not be available to report depending on your type of input sequences or runN vs. runP.
     - |sim_proc_dir| or |sim_overall_dir|
   * - :file:`diamond_annotated_contam.faa/.fnn/.tsv`
     - Generated from EnTAP. Contains all of the transcripts flagged as contaminants (protein/nucleotide format) that are a subset of the diamond_annotated best alignment files. Again this will be seen per database, or overall combining the results from each database used.
     - |sim_proc_dir| or |sim_overall_dir|
   * - :file:`diamond_annotated_without_contam.faa/.fnn/.tsv`
     - Generated from EnTAP. Contains all of the transcripts NOT flagged as contaminants (protein/nucleotide format) that are a subset of the diamond_annotated best alignment files. With this in mind: best_hits = best_hits_no_contam + best_hits_contam. Again this will be seen per database, or overall combining the results from each database used. These sequences are separated from the rest for convenience if you would like to examine them differently
     - |sim_proc_dir| or |sim_overall_dir|
   * - :file:`unannotated.faa/.fnn/.tsv`
     - Generated from EnTAP. Sequences (protein/nucleotide) from the transcriptome that did not hit against this particular reference database. This does not include sequences that were lost during expression filtering or frame selection. Again this will be seen per database, or overall combining the results from each database used.
     - |sim_proc_dir| or |sim_overall_dir|
   * - :file:`diamond_unselected_hits.faa/.fnn/.tsv`
     - Generated from EnTAP. Similarity searching can result in several hits for each query sequence. With only one best alignment being selected, the rest are unselected and end up here. Unselected hits can be due to a low e-value, coverage, or other properties EnTAP takes into account when selecting hits
     - |sim_proc_dir| or |sim_overall_dir|
   * - :file:`species_bar.txt/png`
     - Generated from EnTAP. Bar graph representing the top 10 species that were hit within a reference database or overall. 
       Example:   
	   
       .. image:: plot_sim_species_bar.png
          :scale: 50% 
          :align: center
		  
     - |sim_fig_dir|
   * - :file:`contam_bar.txt/png`
     - Generated from EnTAP. Bar graph representing the top 10 contaminants (within best hits) that were hit against the database or overall.
       Example:   
	   
       .. image:: plot_sim_contam_bar.png
          :scale: 50% 
          :align: center
		  
     - |sim_fig_dir|


Similarity Search Headers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TSV files generated from EnTAP will have the following headers for Similarity Searching (from left to right). Other headers may be present from previous stages of EnTAP (such as Frame Selection or Expression Filtering).

    * Query sequence ID
    * Subject sequence ID
    * Percentage of identical matches
    * Alignment length
    * Number of mismatches
    * Number of gap openings
    * Start of alignment in query
    * End of alignment in query
    * Start of alignment in subject
    * End of alignment in subject
    * Expect (e) value
    * Query coverage
    * Subject title
    * Species (pulled from hit)
    * Origin Database
    * Contaminant (yes/no the hit was flagged as a contaminant)
    * Informative (yes/no the hit is informative)
	
If you ran Similarity Searching against a UniProt database, EnTAP will pull additional UniProt information for your alignments. The following headers will then be added.

    * UniProt Database Cross Reference
    * UniProt Additional Information
    * UniProt KEGG Terms
    * UniProt GO Biological
    * UniProt GO Cellular
    * UniProt GO Molecular