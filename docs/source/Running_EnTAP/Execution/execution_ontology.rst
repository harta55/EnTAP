.. |egg_dir| replace:: :file:`/gene_family/EggNOG`
.. |egg_proc_dir| replace:: :file:`/gene_family/EggNOG/processed`
.. |inter_dir| replace:: :file:`/gene_family/InterProScan`
.. |inter_proc_dir| replace:: :file:`/gene_family/InterProScan/processed`
.. |eggnog_mapper_git| replace:: https://github.com/jhcepas/eggnog-mapper
.. |eggnog_website| replace:: http://eggnog5.embl.de/#/app/home
.. |interproscan_website| replace:: https://www.ebi.ac.uk/interpro/

Gene Family / Ontology Analysis
====================================
This stage of EnTAP attempts to functionally annotate our filtered transcriptome (after Frame Selection and Expression Analysis has been performed) to align information such as gene families, orthologous groups, protein domains, and Gene Ontology terms. EnTAP allows for either :ref:`EggNOG<eggnog-label>` (default and recommended) or :ref:`InterProScan<interproscan-label>` to be used with the :file:`ontology` flag. 

Modify the :file:`ontology` flag within the |run_ini_file_format| by using 0 (EggNOG) and/or 1 (InterProScan) in a comma-separated list.

.. code-block:: bash

    ontology=0
	
To include InterProScan, simply do the following:

.. code-block:: bash

    ontology=0,1

General Ontology Analysis Flags
------------------------------------

.. list-table:: **Ontology Flags**
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1    
   
   * - param
     - description
     - location (cmd/R-ini,E-ini)
     - qualifier
     - example
   * - ontology_source
     - Specify which ontology source packages you would like to use. Multiple flags may be used to specify execution of multiple software packages.
            * 0 - EggNOG (default)
            * 1 - InterProScan
     - R-ini
     - multi-integer
     - 0

.. _eggnog-label:

EggNOG Analysis
-----------------------
By default, EnTAP will utilize EggNOG-mapper (|eggnog_mapper_git|) to access the collection of EggNOG databases (|eggnog_website|) to utilize orthology relationships to assign a myriad of functional information. This is a very powerful tool, especially for non-model transcriptomes where functional data may be limited. 

EggNOG analysis is executed by default with EnTAP so the only thing to make sure of is that the database and execution paths are correct within both ini files. Optional contaminant analysis can be turned on/off, outlined :ref:`here<eggnog_contam-label>`.

.. _eggnog_contam-label:

EggNOG Contaminant Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
EggNOG contaminant analysis can be turned on/off through the :file:`eggnog-contaminant` flag (on by default). When turning EggNOG contaminant analysis on, be sure to review the Similarity Search :file:`contam` flag as the same taxons specified there will be used (must exist in the NCBI Taxonomy Database). 

If EggNOG contaminant analysis is turned on, the results will be displayed in the Log File under the EggNOG section. A contaminant is determined by taking the narrowest Orthologous Groups (seen as "EggNOG Member OGs" in the EnTAP output) assigned to each query and comparing its full lineage to the contaminants input by the user. Although EnTAP reports this information, it does not mean that the query is automatically considered a contaminant through EggNOG analysis. The final contaminant status (seen as "Contaminant" in the EnTAP output) of a query will first be determined through Similarity Search then, if no alignment is found through Similarity Search, rely on the EggNOG contaminant analysis. 

EggNOG Commands
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Ontology - EggNOG Specific Flags
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1    
   
   * - param
     - description
     - location (cmd/R-ini,E-ini)
     - qualifier
     - example
   * - eggnog-map-data
     - Path to the directory containing the EggNOG SQL database |eggnog_map_sql_db_file_format| that was downloaded during the Configuration stage. EnTAP will check for the eggnog.db database within this specified directory
     - E-ini
     - string
     - /path/to/eggnog_db_directory
   * - eggnog-map-dmnd
     - Path to the EggNOG DIAMOND configured database |eggnog_map_dmnd_db_file_format| that was generated during the Configuration stage. 
     - E-ini
     - string
     - /databases/eggnog_proteins.dmnd
   * - eggnog-map-exe
     - Path to the EggNOG-mapper executable, or method of execution. If installed globally, this is simply |emapper_exe_format|
     - E-ini
     - string
     - emapper.py
   * - eggnog-contaminant
     - Specify this to turn on/off EggNOG contaminant analysis (on by default). This leverages the taxon input from the contaminant Similarity Search command to  determine if an EggNOG annotation should be flagged as a contaminant. EggNOG contaminant analysis can only be performed alongside Similarity  Search contaminant analysis (not on its own) and will only be utilized if no alignments were found for a given transcript during Similarity Searching
     - R-ini
     - bool
     - true
   * - eggnog-dbmem
     - Specify this to use the '--dbmem' flag with EggNOG-mapper. This will load the entire eggnog.db sqlite3 database into memory which can require up to ~44GB of memory. However, this will significantly speed up EggNOG annotations
     - R-ini
     - bool
     - true
   * - eggnog-sensitivity
     - Specify the DIAMOND sensitivity used during EggNOG mapper execution against the EggNOG database. Sensitivities are based off of DIAMOND documentation with a higher sensitivity generally taking longer but giving a higher alignment rate. Sensitivity options are fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensitive.
     - R-ini
     - string
     - more-sensitive
	 
Interpreting EggNOG Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The |egg_dir| directory will contain all of the relevant information for the EggNOG stage of the pipeline. This folder will contain the files generated from EggNOG-mapper alongside the files generated by EnTAP. EnTAP files can be found within the |egg_proc_dir| directory.

Below are example files with a transcriptome labelled 'transcriptome' utilizing runP. 

.. list-table:: **EggNOG Results**
   :align: left
   :widths: 10 50 10
   :header-rows: 1    
   
   * - filename
     - description
     - directory
   * - :file:`blastp_transcriptome.emapper.annotations`
     - Generated from EggNOG-mapper. Contains important functional annotation information pulled from orthologous group alignment within EggNOG databases. This will be prepended with blastp or blastx depending on if runP or runN were used.
     - |egg_dir|
   * - :file:`blastp_transcriptome.emapper.seed_orthologs`
     - Generated from EggNOG-mapper. Contains all assigned seed orthologs for the sequences that were ran using EggNOG-mapper. Information in this is similar to that seen with DIAMOND or BLAST runs such as e-value and coverages. This will be prepended with blastp or blastx depending on if runP or runN were used.
     - |egg_dir|
   * - :file:`blastp_transcriptome.emapper.hits`
     - Generated from EggNOG-mapper. Contains all of the hits against the EggNOG database (from DIAMOND). EggNOG-mapper will first align our input transcriptome to the EggNOG database which can result in multiple hits. The selected hits are seen in the .emapper.seed_orthologs file while the rest remain here. This will be prepended with blastp or blastx depending on if runP or runN were used.
     - |egg_dir|
   * - :file:`eggnog_unannotated.fnn/faa`
     - Generated from EnTAP. Sequences where NO alignnment was made with the EggNOG database (nucleotide/protein).
     - |egg_proc_dir|
   * - :file:`eggnog_annotated.fnn/faa`
     - Generated from EnTAP. Sequences where an alignnment was made with the EggNOG database (nucleotide/protein).
     - |egg_proc_dir|
   * - :file:`eggnog_contaminants.fnn/faa`
     - Generated from EnTAP. Sequences that were flagged as a contaminant after EggNOG analysis
     - |egg_proc_dir|

EggNOG Headers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TSV files generated from EnTAP will have the following headers from EggNOG analysis.

    * EggNOG Seed Ortholog
    * EggNOG Seed E-Value
    * EggNOG Seed Score
    * EggNOG Tax Scope Max
    * EggNOG Member OGs
    * EggNOG Description
    * EggNOG COG Abbreviation
    * EggNOG COG Description
    * EggNOG BIGG Reaction
    * EggNOG KEGG KO
    * EggNOG KEGG Pathway
    * EggNOG KEGG Module
    * EggNOG KEGG Reaction
    * EggNOG KEGG RClass
    * EggNOG BRITE
    * EggNOG GO Biological
    * EggNOG GO Molecular
    * EggNOG Protein Domains
    * Contaminant

.. _interproscan-label:

InterProScan Analysis
-------------------------
The user has the option to use InterProScan (|interproscan_website|) as an additional method of determining functional annotation of our transcripts. InterProScan is a powerful tool that will classify our transcripts into families to predict domains and other important functional information.

Running InterProScan
^^^^^^^^^^^^^^^^^^^^^^^^
In order to run InterProScan, as mentioned above, the :file:`ontology` flag must also include '1' within the |run_ini_file_format| file. Additional parameters can be set, such as which additional database should be analyzed through the :file:`protein` command. These can be seen below.

InterProScan Commands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. list-table:: Ontology - InterProScan Specific Flags
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1    
   
   * - param
     - description
     - location (cmd/R-ini,E-ini)
     - qualifier
     - example
   * - interproscan-db
     - User this option if you would like to run InterProScan against specific databases. Multiple databases can be selected. 
           * tigrfam
           * sfld
           * prodom
           * hamap
           * pfam
           * smart
           * cdd
           * prositeprofiles
           * prositepatterns
           * superfamily
           * prints
           * panther
           * gene3d
           * pirsf
           * coils
           * mobidblite
     - R-ini
     - multi-string
     - pfam
   * - interproscan-exe
     - Specify the execution method for InterProScan. Commonly this can be the path to the :file:`interproscan.sh` file
     - E-ini
     - string
     - interproscan.sh
	 
	
Interpreting InterProScan Results
-------------------------------------------
The |inter_dir| directory will contain all of the relevant information for the optional InterProScan stage of the pipeline. This folder will contain files generated by InterProScan as well as those by EnTAP (|inter_proc_dir|).

Below are the example files you may find when including InterProScan:

.. list-table:: **InterProScan Results**
   :align: left
   :widths: 10 50 10
   :header-rows: 1    
   
   * - filename
     - description
     - directory
   * - :file:`interproscan.tsv/xml`
     - Generated from InterProScan. Tab delimited or XML file containing information on the sequences with domain matches. Information such as signature accession/description information and GO/Pathway alignments.
     - |inter_dir|
   * - :file:`unannotated_sequences.fnn/faa`
     - Generated from EnTAP. Sequences where NO domain could be assigned (nucleotide/protein) through InterProScan
     - |inter_proc_dir|
   * - :file:`annotated_sequences.fnn/faa`
     - Generated from EnTAP. Sequences where a domain could be assigned (nucleotide/protein) through InterProScan
     - |inter_proc_dir|

InterProScan Headers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TSV files generated from EnTAP will have the following headers from InterProScan analysis.

    * IPScan GO Biological
    * IPScan GO Cellular
    * IPScan GO Molecular
    * IPScan Pathways
    * IPScan InterPro ID
    * IPScan Protein Database
    * IPScan Protein Description
    * IPScan E-Value
