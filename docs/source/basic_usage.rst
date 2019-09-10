.. _NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/taxonomy
.. _Bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. |out_dir| replace:: /entap_outfiles
.. |libs_dir| replace:: /libs
.. |entap_dir| replace:: /EnTAP
.. |src_dir| replace:: /src
.. |config_file| replace:: entap_config.txt
.. |bin_dir| replace:: /bin
.. |test_dir| replace:: /test_data
.. |data_dir| replace:: /databases
.. |tax_file| replace:: download_tax.pl
.. |graph_file| replace:: entap_graphing.py
.. |go_term| replace:: go_term.entp
.. |tax_bin| replace:: ncbi_tax_bin.entp
.. |tax_data| replace:: ncbi_tax.entp

.. |ref_comp| replace:: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/
.. |ref_plant| replace:: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/
.. |ref_mamm| replace:: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/
.. |ref_nr| replace:: ftp://ftp.ncbi.nlm.nih.gov/blast/db/
.. |uni_swiss| replace:: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
.. |uni_trembl| replace:: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
.. |entap_bin_ftp| replace:: https://treegenesdb.org/FTP/EnTAP/latest/databases/entap_database.bin.gz
.. |entap_sql_ftp| replace:: https://treegenesdb.org/FTP/EnTAP/latest/databases/entap_database.db.gz
.. |eggnog_sql_ftp| replace:: http://eggnogdb.embl.de/download/emapperdb-4.5.0/eggnog.db.gz
.. |eggnog_fasta_ftp| replace:: http://eggnogdb.embl.de/download/latest/eggnog-mapper-data/eggnog4.clustered_proteins.fa.gz

.. |flag_path| replace:: paths
.. |flag_taxon| replace:: taxon


Basic Usage
============

*EnTAP* has two stages, :ref:`configuration<config-label>` and :ref:`execution<exe-label>`. Configuration will download and configure databases that are required for execution. It  will need to be completed before the first execution and if any of the source databases would like to be updated. Execution is the main stage of EnTAP that will annotate a transcriptome input by the user. 

.. _config-label:

Configuration
-------------
Configuration is the first stage of EnTAP that will download and configure the necessary databases for full functionality. This is run if you would like to change/update the databases that EnTAP is reading from. I'll break this up into two sections, :ref:`Config File<hierarchy-label>` and :ref:`Usage<usage-label>`. The Config File section will just describe how to ensure EnTAP is reading from the correct paths, which can be easily changed in the |config_file| (more on that later!). It will also go over the directories included in the installation. The Usage sections will go over the basic usage during the Configuration stage of EnTAP and how to setup reference databases. 


.. _hierarchy-label:

Config File
^^^^^^^^^^^^^^^^^

From here on out, the "execution", or "EnTAP", directory will refer to the directory containing the EnTAP install (or binary file). Typically, this will just be at the root directory that was downloaded from the repository. All paths mentioned in this documentation will be relative to this directory. 


Why is this important? EnTAP relies on several accompanying software packages and databases in order to run properly. Correct recognition of these paths is crucial and, as such, needed an entire section! The |config_file| is the answer to this pathing issue. It contains all of the necessary paths required for EnTAP to run and can be configured as seen fit. 

Normally when a user is trying to execute EnTAP, they will need to specify the path to this Config file with the |flag_path| flag. It's never a bad idea to specify this. However, in an attempt to make things a bit easier for the end user, EnTAP can "assume" some default paths so --|flag_path| does not always need to be specified. If the flag is not specified , EnTAP will check if the Config File resides in the current working directory (not execution directory). If an empty |config_file| file is found in the working directory, the following assumptions will be made as far as paths are concerned:

    * diamond_exe_path=/EnTAP/libs/diamond-0.9.9/bin/diamond
    * rsem_exe_path=/EnTAP/libs/RSEM-1.3.0 (this is a path to the directory)
    * genemarkst_exe_path=/EnTAP/libs/gmst_linux_64/gmst.pl
    * eggnog_sql_database=/EnTAP/databases/eggnog.db (downloaded during Configuration)
    * eggnog_dmnd_database=/EnTAP/bin/eggnog_proteins.dmnd (downloaded during Configuration)
    * interpro_exe_path=interproscan.sh
    * entap_database_sql_path=/EnTAP/databases/entap_database.db (text version, downloaded during Configuration)
    * entap_database_bin_path=/EnTAP/bin/entap_database.bin (binary version, downloaded during Configuration)
    * entap_graphing_script=/EnTAP/src/entap_graphing.py


Keep in mind, it is always safer to use the |flag_path| flag to specify the |config_file|. EnTAP will always recognize these paths first and they will override defaults above. 

If something is globally installed, such as "interpro_exe_path" above, put how you'd normally run the software after the '='. As another example, running DIAMOND through a global installation may simply be "diamond". The Config File line for DIAMOND will simply read:

    * diamond_exe_path=diamond

.. note:: Either the EnTAP binary database (default) or the EnTAP SQL database is required for execution. Both are not needed.

.. warning:: Be sure you set the paths before moving on (besides the databases that haven't been downloaded yet)! 


.. _usage-label:

Usage - Preparing Your Reference Databases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All source databases must be provided in FASTA format (protein) so that they can be indexed for use by DIAMOND.  This can be completed independent of EnTAP with DIAMOND (- - makedb flag) or as part of the Configuration phase of EnTAP. This section will focus on downloading and preparing some of the more common FASTA source databases. If you already have DIAMOND databases configured, you can skip to :ref:`Usage - Running Configuration<usage_config-label>`. Even if you have a DIAMOND database already configured, Configuration must still be ran!


While any protein FASTA database can be used, it is recommended to use NCBI (Genbank) sourced databases such as RefSeq databases or NR.  In addition, EnTAP can easily accept EBI databases such as UniProt/SwissProt.  

EnTAP can recognize the species information from these header formats ONLY (NCBI and UniProt):

* [homo sapiens]

* OS=homo sapiens

If the individual FASTAs in a custom database you create do not adhere to one of these two formats, it will not be possible to weight taxonomic or contaminant status from them. You will need to change the headers to ensure they align. 

The following FTP sites contain common reference databases that EnTAP can recognize:
   * RefSeq: |ref_comp|

   * Plant RefSeq: |ref_plant|

   * Mammalian RefSeq: |ref_mamm|

   * NR: |ref_nr|

   * SwissProt: |uni_swiss|
   
       * Reviewed
       * It is highly recommended to use the UniProt SwissProt database as EnTAP will map all UniProt alignments to additional database cross-references

   * TrEMBL: |uni_trembl|
   
       * Unreviewed

Both Uniprot databases (SwissProt and TrEMBL) can be downloaded on a Unix system through the following command:

.. code-block:: bash
 
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

Or, for the TrEMBL database:

.. code-block:: bash

    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

Alternatively, the NCBI databases must be downloaded in separate, smaller files, and concatenated together. As an example, the following commands will download and combine the NR database files:

Download:

.. code-block:: bash

    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz

Decompress/Concatenate:

.. code-block:: bash

    tar -xvzf nr.*
   
    cat nr.* > nr_database.fasta
    

....

It is generally recommended that a user select at least three databases with varying levels of curation.  Unless the species is very non-model (i.e. does not have close relatives in databases such as RefSeq, it is not necessary to use the full NR database which is less curated). Once your FASTA databases are ready, move on to :ref:`Running Configuration<usage_config-label>`.


.. _usage_config-label:

Usage - Running Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have your protein FASTA database ready, you can begin to run the Configuration stage. As mentioned before, Configuration will only need to be run once prior to :ref:`Execution<exe-label>` unless you would like to configure/update more databases. 

To run configuration with a FASTA database to output directory path/to/output (default is current working directory), the command is as follows (additional databases can be specified if necessary with the -d flag and threads with the -t flag):

.. code-block:: bash

    EnTAP --config -d path/to/database.fasta -d path/to/database2.fasta --out-dir path/to/output -t 8


Configuration can be run without formatting a FASTA database for DIAMOND is as follows with 8 threads:

.. code-block:: bash

    EnTAP --config -t 8

.. note:: This is the only stage that requires connection to the Internet.

In both cases, the following databases will be downloaded and configured:

* EnTAP Binary Database:
    * Comprised of Gene Ontology, UniProt, and Taxonomic mappings for use during Execution. FTP downloaded file.
    * Downloaded from |entap_bin_ftp|
    * Filename: entap_database.bin
    * The SQL version is the same database, but formatted as a SQL database. Only one version of the database is needed (binary is used by default)

* EggNOG DIAMOND Reference:
    * Reference database containing EggNOG database entries
    * FASTA file is downloaded and configured for DIAMOND from |eggnog_fasta_ftp|
    * Filename: eggnog_proteins.dmnd

* EggNOG SQL Database:
    * SQL database containing EggNOG mappings
    * Downloaded from |eggnog_sql_ftp|
    * Filename: eggnog.db

The EnTAP Binary Database is downloaded from the FTP addresses below. By default, the binary version will be downloaded and used. Only one version is required. If you experience any trouble in downloading, you can simply specify the - - data-generate flag during Configuration to configure it locally (more on that later). The database for the newest version of EnTAP will always reside in the "latest" FTP directory. Keep in mind, if you are using an older version of EnTAP, you do not want to download from the "latest" directory. Instead, you will need to consider the version you are using. The FTP will always be updated only when a new database version is created. For example, if you see v0.8.2 and v0.8.5 on the FTP while you are using v0.8.3, you will download the database located in the v0.8.2 directory. 

    * |entap_bin_ftp|
    * |entap_sql_ftp|


.. warning ::
    DIAMOND databases must be configured and eventually executed with the same version of DIAMOND.

Configuration Flags:
^^^^^^^^^^^^^^^^^^^^^^

Required Flags:

* (- - config)
    * The only required flag. 
    * Although in order to run the full EnTAP pipeline, you must have a .dmnd configured database.

Optional Flags:

* (-d/ - - database)
    * Specify any number of FASTA formatted databases you would like to configure for EnTAP
    * Not necessary if you already have DIAMOND configured databases (.dmnd)

* (- - |flag_path|)
    * Point to |config_file| to specify file paths
    * DIAMOND is the only path necessary during Configuration
    * Default: |config_file| residing in the current working directory

* (- -  out-dir)
    * Specify an output directory for the databases to be sent to (recommended)
    * This will send the EnTAP database and DIAMOND databases to this location

* (- t/ - - threads)
    * Specify thread number for Configuration

* (- - data-generate)
    * Specify this flag is you would like to generate the EnTAP database rather than downloading from FTP (default)
    * I'd only use this if you're having issues with the FTP

* (- - data-type)
    * Specify which databases you'd like to generate/download

        * 0. Binary Database (default) - This will be much quicker and is recommended
        * 1. SQL Database - Slower although will be more easily compatible with every system

    * This can be flagged multiple times (ex: - - data-type 0 - - data-type 1)
    * I would not use this flag unless you are experiencing issues with the EnTAP Binary Database

.. test-label:

Test Data
-------------
Before continuing on to the :ref:`Execution<exe-label>` stage, it is advised to do a test run of EnTAP to ensure that everything is properly configured. There should be no errors in the test run. The test data resides within the |test_dir| directory of the main EnTAP directory. This will walk you through configuring a database for DIAMOND (if you haven't already done so) and executing EnTAP with and without frame selection. 

Before we begin, make sure that the paths in the configuration file are correct. Since we are running the configuration stage, EnTAP will check to make sure you have the other databases downloaded (which should have been done prior to this). To begin the test, execute the following command to configure the test DIAMOND database:

.. code-block:: bash

    EnTAP --config -d /test_data/swiss_prot_test.fasta --out-dir /test_data


This should finish very shortly without any errors and you should find a swiss_prot_test.dmnd file within the |test_dir| directory. 

Next up is verifying the main execution stage! Once again, first ensure that the Config File has all of the correct paths. We are going to check an execution with and without frame selection. If you are not going to use frame selection, you may skip this test!

.. note:: The following tests will take longer as they will be testing the entire pipeline and running against the larger EggNOG database.

To test EnTAP with frame selection, execute the following command:

.. code-block:: bash

    EnTAP --runP -i /test_data/trinity.fnn -d /test_data/bin/swiss_prot_test.dmnd

To test EnTAP without frame selection, execute the following command:

.. code-block:: bash

    EnTAP --runP -i /test_data/trinity.faa -d /test_data/swiss_prot_test.dmnd

These should run without error and you should have several files within the created |out_dir| directory. The final_annotations_lvl0.tsv file should resemble the test_data/final_annotations_test.tsv file. 

If any failures were seen during the above executions, be sure to go through each stage of installation and configuration to be sure everything was configured correctly before continuing!

.. _exe-label:

Execution
-------------
The Execution stage of EnTAP is the main annotation pipeline. After Configuration is run at least once, this can be run continually without requiring Configuration to be ran again (unless more databases will be configured). 

The following stages will be run:

#. :ref:`Expression Filtering<exp-label>` (optional)
#. :ref:`Frame Selection<frame-label>` (optional)
#. Similarity Search
#. Orthologous Group Assignment
#. InterProScan (optional)

Input Files:
^^^^^^^^^^^^^^^^^
Required:

* .FASTA formatted transcriptome file (either protein or nucleotide)
* .dmnd (DIAMOND) indexed databases, which can be formatted in the :ref:`Configuration<config-label>` stage. 

Optional:

* .BAM/.SAM alignment file. If left unspecified expression filtering will not be performed. 
    * This can be generated by software that does not perform gapped alignments such as `Bowtie`_ (not Bowtie2). All you need to generate an alignment file is a pair of reads and your assembled transcriptome!

Sample Run:
^^^^^^^^^^^^^^^^^

A specific run flag (**runP/runN**) must be used:

* runP: Indicates blastp. Frame selection will be ran if nucleotide sequences are input
* runN: Indicates blastx. Frame selection will not be ran with this input


An example run with a nucleotide transcriptome (transcriptome.fasta), two reference DIAMOND databases, an alignment file (alignment.sam), and 8 threads:

.. code-block:: bash

    EnTAP --runP -i path/to/transcriptome.fasta -d path/to/database.dmnd -d path/to/database2.dmnd -a path/to/alignment.sam -t 8


With the above command, the entire EnTAP pipeline will run. Both Frame Selection and Expression Filtering can be skipped if preferred by the user. If a protein transcriptome is input with the runP flag, or the runN flag is used, Frame Selection will be skipped.  If there is not a short read alignment file provided in SAM/BAM format, then Expression Filtering via RSEM will be skipped. 


Flags:
^^^^^^^^^^^^^^^^^^^^^

Required Flags:

* (- - runP/- - runN)
    * Specify a blastp or blastx annotation
    * If - -runP is selected with a nucleotide input, frame selection will be ran and annotation stages will be executed with protein sequences (blastp)
    * If - -runP is selected with a protein input, frame selection will not be ran and annotation will be executed with protein sequences (blastp)
    * If - -runN is selected with nucleotide input, frame selection will not be ran and annotation will be executed with nucleotide sequences (blastx)

* (-i/- - input)
    * Path to the transcriptome file (either nucleotide or protein)

* (-d/- - database)
    * Specify up to 5 DIAMOND indexed (.dmnd) databases to run similarity search against

Optional Flags:

* (- - |flag_path|)
    * Point to |config_file| to specify proper database and execution paths
    * Default: |config_file| residing in the current working directory

* (-a/- -align)
    * Path to alignment file (either SAM or BAM format)
    * **Note:** Ignoring this flag will skip expression filtering
    * If you have ran alignment with single end reads be sure to use the - -single-end flag as well (paired-end is default)
    * Be sure to specify an FPKM threshold

* (- - contam)
    * Specify :ref:`contaminant<tax-label>` level of filtering
    * Multiple contaminants can be selected through repeated flags

* (- - taxon)
    * This flag will allow for :ref:`taxonomic<tax-label>` 'favoring' of hits that are closer to your target species or lineage. Any lineage can be used as referenced by the NCBI Taxonomic database, such as genus, phylum, or species.
    * Format **must** replace all spaces with underscores ('_') as follows: "- -taxon homo_sapiens" or "- -taxon primates"

* (- - level)
    * Specify Gene Ontology levels you would like to normalize to
    * Any amount of these flags can be used
    * Default: 0 (every level), 3, 4
    * More information at: http://geneontology.org/page/ontology-structure

* (- - out-dir)
    * Specify output folder labelling.
    * Default: /outfiles

* (- - fpkm)
    * Specify FPKM cutoff for expression filtering
    * Default: 0.5

* (-e)
    * Specify minimum E-value cutoff for similarity searching
    * Default: 10E-5

* (- - tcoverage)
    * Specify minimum target coverage for similarity searching
    * Default: 50%

* (- - qcoverage)
    * Specify minimum query coverage for similarity searching
    * Default: 50%

* (- - overwrite)
    * All previously ran files will be overwritten if the same - -tag flag is used
    * Without this flag EnTAP will :ref:`recognize<over-label>` previous runs and skip things that were already ran

* (- - single-end)
    * Signify your reads are single end for RSEM execution
    * Default: paired-end 

* (- - graph)
    * This will check whether or not your system has graphing functionality supported and exit
    * If Python with the Matplotlib module are installed on your system graphing should be enabled!
    * This can be specified on its own

* (-t/ - - threads)
    * Specify the number of threads of execution

* ( - - trim)
    * This flag will trim your sequence headers to anything before a space. It will make your data easier to read if you have a lot of excess information you do not need in your headers.
    * Example: 
   
        * >TRINITY_231.1 protein12312_43_inform
        * >TRINITY_231.1

* (- - state)
    * Precise control over execution :ref:`stages<state-label>`. This flag allows for certain parts to be ran while skipping others. 
    * Warning: This may cause issues depending on what you plan on running! 

* (- - ontology)
    * Specify which ontology packages you would like to use

        * 0 - EggNOG (default)
        * 1 - InterProScan

    * Both or either can be specified with multiple flags

        * Ex: - - ontology 0 - - ontology 1
        * This will run both EggNOG and InterProScan 

* (- - protein)
    * Use this option if you would like to run InterProScan
    * Specify databases to run against (you must have them already installed)
      
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

* (- - version)
    * Prints the current EnTAP version you are running

* (- - uninformative)
    * Path to a list of terms you would like to be deemed "uninformative"
    * The file **must** be formatted with one term on each line of the file
    * Example (defaults):
    
        * conserved
        * predicted
        * unnamed
        * hypothetical
        * putative
        * unidentified
        * uncharacterized
        * unknown
        * uncultured
        * uninformative

* (- - no-check)
    * EnTAP checks execution paths and inputs prior to annotating to prevent finding out your input was wrong until midway through a run. Using this flag will eliminate the check (not advised to use!)

* (- - output-format)
    * Specify multiple output file formats for each stage of the pipeline

        * 1. TSV File (default)
        * 2. CSV File
        * 3. FASTA Protein File (default)
        * 4. FASTA Nucleotide File (default)

* (- - data-type)
    * Specify which database you'd like to execute against (not advised to use)

        * 0. Binary Database (default) - This will be much quicker and is recommended
        * 1. SQL Database - Slower 

    * If you flag this multiple times during execution, EnTAP will just select the first one you input


.. _exp-label:

Expression Analysis
^^^^^^^^^^^^^^^^^^^^^^^
The goal of expression filtering, or transcript quantification, is to determine the relative 
abundance levels of transcripts when taking into account the sequenced reads and how they map 
back to the assembled transcriptome and using this information to filter out suspect expression 
profiles possibly originated from poor or incomplete assemblies. Filtering is done through the use
of the FPKM (fragments per kilobase per of million mapped reads) , or a measurable number of 
expression. This can be specified with the - -fpkm flag as specified above. EnTAP will use this FPKM value
and remove any sequences that are below the threshold.

.. _frame-label:

Frame Selection
^^^^^^^^^^^^^^^^^^
Frame selection is the process of determining the coding region of a transcript. Oftentimes, due to 
assembly errors or other factors, a coding region may not be found for a transcript and EnTAP will remove
this sequence. When a coding region is found, EnTAP will include the sequence for further annotation.

.. _tax-label:

Taxonomic Favoring and Contaminant Filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Taxonomic contaminant filtering (as well as taxonomic favoring) is based upon the `NCBI Taxonomy`_ database. In saying this, all species/genus/lineage names must be contained within this database in order for it to be recognized by EnTAP. 

**Contaminant Filtering:**

Contaminants can be introduced during collection or processing of a sample. A contaminant is essentially a species that is not of the target species you are collecting. Some common contaminants are bacteria and fungi that can sometimes be found within collected samples. If a query sequence from your transcriptome is found when matching against a similarity search database, it will be flagged as such (but NOT removed automatically). Oftentimes, researchers would like to remove these sequences from the dataset. 

An example of flagging bacteria and fungi as contaminants can be seen below:

.. code-block:: bash

    EnTAP --runP -i path/to/transcriptome.fasta -d path/to/database.dmnd -c fungi -c bacteria


**Taxonomic Favoring**

During best hit selection of similarity searched results, taxonomic consideration can utilized. If a certain lineage (such as sapiens) is specified, hits closer in taxonomic lineage to this selection will be chosen. Any lineage such as species/kingdom/phylum can be utilized as long as it is contained within the Taxonomic Database. If it is not located within the database, EnTAP will stop the execution immediately and let you know! 

This feature can be utilized with the |flag_taxon| flag. An example command utilizing both common contaminants and a species taxon can be seen below:

.. code-block:: bash

    EnTAP --runP -i path/to/transcriptome.fasta -d path/to/database.dmnd -c fungi -c bacteria --taxon sapiens

Keep in mind, EnTAP will weigh the E-Value (within a database)and Coverage of the alignment before taxonomic weight in order to provide the most accurate result. If both the E-Value and Coverage are relatively similar, EnTAP will leverage taxonomic information.

.. _over-label:

Picking Up Where You Left Off
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to save time and make it easier to do different analyses of data, EnTAP allows for picking up where you left off if certain stages were already ran and you'd like analyze data with different contaminant flags or taxonomic favoring. As an example, if similarity searching was ran previously you can skip aligning against the database and analyze the data to save time. However, the - - overwrite flag will not allow for this as it will remove previous runs and not recognize them. 

In order to pick up and skip re-running certain stages again, the files that were ran previously **must** be in the same directories and have the same names. With an input transcriptome name of 'transcriptome' and example database of 'complete.protein':

* Expression Filtering
    * transcriptome.genes.results

* Frame Selection
    * transcriptome.fasta.faa
    * transcriptome.fasta.fnn
    * transcriptome.fasta.lst

* Similarity Search
    * blastp_transcriptome_complete.protein.faa.out

* Gene Family
    * blastp_transcriptome_eggnog_proteins.out (for runP)
    * blastp_transcriptome_eggnog_proteins.out (for runN)


Since file naming is based on your input as well, the flags below **must** remain the same:

* (- - runN / - - runP)

* (- - ontology)

* (- - protein)

* (-i / - - input)

* (-a / - - align)

* (-d / - - database)
    * Does not necessarily need to remain the same. If additional databases are added, EnTAP will recognize the new ones and run similarity searching on them whilst skipping those that have already been ran

* (- - qcoverage)

* (- - tcoverage)

* (- - trim)

* (- - out-dir)


.. _state-label:

State Control
^^^^^^^^^^^^^^

.. warning:: This is experimental and certain configurations may not work. This is not needed if you'd like to run certain portions because of "picking up where you left off!"

State control of EnTAP allows you to further customize your runs. This is separate from the exclusion of - - align flag to skip expression filtering, or runP, instead of runN, to skip frame selection. You probably will never actually have to use this feature! Nonetheless, state control is based around the following stages of EnTAP:

#. Expression Filtering
#. Frame Selection
#. Transcriptome Filtering (selection of final transcriptome)
#. Similarity Search
#. Gene Ontology / Gene Families

With this functionality of EnTAP, you can execute whatever states you would like with certain commands. Using a '+' will execute from that state to the end, while using a 'x' will stop at that state. These basic commands can be combined to execute whatever you would like. It's easier if I lay out some examples:

* (- - state 1+)
    * This will start at expression filtering and continue to the end of the pipeline

* (- - state 1+4x)
    * This will start at expression filtering and stop after similarity search

* (- - state 4x)
    * This will just execute similarity search and stop

* (- - state 1+3x5)
    * This will essentially execute every stage besides similarity searching

The default 'state' of EnTAP is merely '+'. This executes every stage of the pipeline (or attempts to if the correct commands are in place). 
