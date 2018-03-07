.. _NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/taxonomy
.. _Bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. |out_dir| replace:: /outfiles
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

.. |flag_path| replace:: paths
.. |flag_taxon| replace:: taxon


Basic Usage
============

*EnTAP* has two stages of execution, :ref:`configuration<config-label>` and :ref:`run<run-label>`. Configuration should be completed before the first run and everytime any of the source databases have been updated by the user.  This should also be run if you would like to include the latest version of the NCBI Taxonomy database and the Gene Ontology database.  All of these are updated regularly at the source and you can ensure you have the most recent version by running configuration before your annotation runs.

.. _config-label:

Configuration
-------------
Configuration is the first stage of EnTAP that will download and configure the necessary databases for full functionality. This is run if you would like to change/update the databases that EnTAP is reading from. I'll break this up into two sections, :ref:`folder hierarchy<hierarchy-label>` and :ref:`usage<usage-label>`. The folder hierarchy section will just describe how everything is 'typically' setup with EnTAP, however these paths can be easily changed in the |config_file| (more on that later!). The usage section will go over the basic usage during the Configuration stage of EnTAP. 


.. _hierarchy-label:

Folder Hierarchy
^^^^^^^^^^^^^^^^^

The EnTAP folder organization is referred to as the execution directory where all files will be made available. This is essentially the hierarchy that was downloaded from the repository. 

The |entap_dir| directory contains:

    * |entap_dir| |libs_dir| 
    * |entap_dir| |src_dir|
    * |entap_dir| |bin_dir| (created during configuration)
    * |entap_dir| |data_dir| (created during configuration)

In addition to some other files/directories.

Recognition of EnTAP databases, src files, and execution accompanying pipeline software can rely on this 'default' directory hierarchy. However, any necessary files/directories can be changed from the default with the  |config_file| file (by specifying the |flag_path| flag). 

.. warning:: Ensure you are pointing to the correct paths if not using the defaults!

The |config_file| file mentioned above has the following defaults:

    * diamond_exe_path=/EnTAP/libs/diamond-0.8.31/bin/diamond
    * rsem_exe_path=/EnTAP/libs/RSEM-1.3.0 (this is a path to the directory)
    * genemarkst_exe_path=/EnTAP/libs/gmst_linux_64/gmst.pl
    * eggnog_exe_path=/EnTAP/libs/eggnog-mapper/emapper.py
    * eggnog_download_exe=/EnTAP/libs/eggnog-mapper/download_eggnog_data.py
    * eggnog_database=/EnTAP/libs/eggnog-mapper/data/eggnog.db (downloaded during Configuration)
    * entap_tax_database=/EnTAP/bin/ncbi_tax_bin.entp (binary version, downloaded during Configuration)
    * entap_tax_download_script=/EnTAP/src/download_tax.pl
    * entap_go_database=/EnTAP/bin/go_term.entp (binary version, downloaded during Configuration)
    * entap_graphing_script=/EnTAP/src/entap_graphing.py
    * interpro_exe_path=interproscan.sh


These can be changed to whichever path you would prefer. If something is globally installed, just put how you'd normally run the software after the '=', such as 'diamond' for DIAMOND. EnTAP will recognize these paths first and they will override defaults. 


This configuration file will be automatically detected if it is in the same directory as the EnTAP .exe, otherwise the path to it can be specified through the |flag_path| flag. 

.. note:: Be sure you set the paths before moving on (besides the databases that haven't been downloaded yet)!


.. _usage-label:

Usage
^^^^^

All source databases must be provided in FASTA format so that they can be indexed for use by DIAMOND.  This can be completed independent of EnTAP with DIAMOND (- - makedb flag) or as part of the configuration phase of EnTAP.  While any FASTA database can be used, it is recommended to use NCBI (Genbank) sourced databases such as RefSeq databases or NR.  In addition, EnTAP can easily accept EBI databases such as UniProt/SwissProt.  

EnTAP can recognize the species information from these header formats ONLY (NCBI and UniProt):

* [homo sapiens]

* OS=homo sapiens

If the individual FASTAs in a custom database you create do not adhere to one of these two formats, it will not be possible to weight taxonomic or contaminant status from them.  

The following FTP sites contain common reference databases that EnTAP can recognize:
   * RefSeq: |ref_comp|

   * Plant RefSeq: |ref_plant|

   * Mammalian RefSeq: |ref_mamm|

   * NR: |ref_nr|

   * SwissProt: |uni_swiss|
   
       * Reviewed

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

It is generally recommended that a user select at least three databases with varying levels of curation.  Unless the species is very non-model (i.e. does not have close relatives in databases such as RefSeq, it is not necessary to use the full NR database which is less curated).


To run configuration with a sample database, the command is as follows:

.. code-block:: bash

    EnTAP --config -d path/to/database

This stage must be done at least once prior to :ref:`running<run-label>`. Once the database is configured, you need not do it again unless you updated your original database or plan on configuring several others.


.. note:: If you already have DIAMOND (.dmnd) configured databases, you can skip the configuration of that database. Although, due to other EnTAP database downloading (taxonomy and ontology), configuration must still be ran at least once without any flags.

Configuration can be ran without formatting a database as follows:

.. code-block:: bash

    EnTAP --config

In both cases, the following databases will be downloaded:

* NCBI Taxonomic Database (indexed for EnTAP)
* Gene Ontology Database (indexed for EnTAP)
* EggNOG DIAMOND Database
* EggNOG SQL Database

.. note:: This is the only stage that requires connection to the Internet.

If you experience any trouble in downloading the databases indexed for EnTAP (taxonomy and gene ontology), you can use the databases contained in the repo download, databases.tar.gz. Just be sure to set the configuration file to these database paths (as these are the binaries)!

Flags:
^^^^^^^^^^^^^^^^^^^^^

Required Flags:

* (- - config)
    * The only required flag. 
    * Although in order to run the full EnTAP pipeline, you must have a .dmnd configured database.

Optional Flags:

* (-d/ - - database)
    * Specify any number of FASTA formatted databases you would like to configure for EnTAP
    * Not necessary if you already have DIAMOND configured databases (.dmnd)

* (- - |flag_path|)
    * Point to |config_file| for specifying paths

* (- - database-out)
    * Specify an output directory for the databases to be sent to
    * This will send the Taxonomic Database, GO Database, and any DIAMOND databases to this location
    * EggNOG database will not be sent here as it must remain in the EggNOG directory

* (- t/ - - threads)
    * Specify thread number for Configuration

.. test-label:

Test Data
-------------
Before continuing on to the :ref:`run<run-label>` stage, it is advised to do a test run of EnTAP to ensure that everything is properly configured. There should be no errors in the test run. The test data resides within the |test_dir| directory of the main EnTAP directory. This will walk you through configuring a database for DIAMOND (if you haven't already done so) and executing EnTAP with and without frame selection. 

Before we begin, make sure that the paths in the configuration file are correct. Since we are running the configuration stage, EnTAP will check to make sure you have the other databases downloaded (which should have been done prior to this). To begin the test, execute the following command to configure the test DIAMOND database:

.. code-block:: bash

    EnTAP --config -d /test_data/swiss_prot_test.fasta --database-out /test_data


This should finish very shortly without any errors and you should find a uniprot_sprot_test.dmnd file within the |test_dir| directory. 

Next up is verifying the main execution stage! Once again, first ensure that the configuration file has all of the correct paths. We are going to check an execution with and without frame selection. If you are not going to use frame selection, you may skip this test!

.. note:: The following tests will take longer as they will be testing the entire pipeline and running against the larger EggNOG database.

To test EnTAP with frame selection, execute the following command:

.. code-block:: bash

    EnTAP --runP -i /test_data/trinity.fnn -d /test_data/uniprot_sprot_test.dmnd

To test EnTAP without frame selection, execute the following command:

.. code-block:: bash

    EnTAP --runP -i /test_data/trinity.faa -d /test_data/uniprot_sprot_test.dmnd

These should run without error and you should have several files within the created |out_dir| directory. The final_annotations_lvl0.tsv file should resemble the test_data/final_annotations_test.tsv file. 

If any failures were seen during the above executions, be sure to go through each stage of installation and configuration to be sure everything was configured correctly before continuing!

.. _run-label:

Run
-------------
The run stage of *EnTAP* is the main annotation pipeline. After configuration is ran at least once, this can be ran continually without requiring configuration to be ran again (unless more databases will be configured). 

The following stages will be ran:

#. :ref:`Expression Filtering<exp-label>` (optional)
#. :ref:`Frame Selection<frame-label>` (optional)
#. Similarity Search
#. Orthologous Group Assignment
#. InterProScan (optional)

Input Files:
^^^^^^^^^^^^
Required:

* .FASTA formatted transcriptome file (either protein or nucleotide)
* .dmnd (DIAMOND) indexed databases, which can be formatted in the :ref:`configuration<config-label>`stage. 

Optional:

* .BAM/.SAM alignment file. If left unspecified expression filtering will not be performed. 
    * This can be generated by software that does not perform gapped alignments such as `Bowtie`_ (not Bowtie2). All you need to generate an alignment file is a pair of reads and your assembled transcriptome!

Sample Run:
^^^^^^^^^^^

A specific run flag (**runP/runN**) must be used:

* runP: Indicates blastp. Frame selection will be ran if nucleotide sequences are inputted
* runN: Indicates blastx. Frame selection will not be ran with this input


An example run with a nucleotide transcriptome:

.. code-block:: bash

    EnTAP --runN -i path/to/transcriptome.fasta -d path/to/database.dmnd -d path/to/database2.dmnd -a path/to/alignment.sam


With the above command, the entire EnTAP pipeline will run. Both frame selection and expression filtering can be skipped if preferred by the user.  EnTAP would require protein sequences (indicated by --runP) in order to avoid frame selection.  If there is not a short read alignment file provided in SAM/BAM format, then expression filtering via RSEM will be skipped. 


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
    * This will check whether or not your system has graphing functionality supported
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Taxonomic contaminant filtering (as well as taxonomic favoring) is based upon the `NCBI Taxonomy`_ database. In saying this, all species/genus/lineage names must be contained within this database in order for it to be recognized by EnTAP. 

**Contaminant Filtering:**

Contaminants can be introduced during collection or processing of a sample. A contaminant is essentially a species that is not of the target species you are collecting. Some common contaminants are bacteria and fungi that can sometimes be found within collected samples. If a query sequence from your transcriptome is found when matching against a similarity search database, it will be flagged as such (but NOT removed automatically). Oftentimes, researchers would like to remove these sequences from the dataset. 

An example of flagging bacteria and fungi as contaminants can be seen below:

.. code-block:: bash

    EnTAP --runN -i path/to/transcriptome.fasta -d path/to/database.dmnd -c fungi -c bacteria


**Taxonomic Favoring**

During best hit selection of similarity searched results, taxonomic consideration can utilized. If a certain lineage (such as sapiens) is specified, hits closer in taxonomic lineage to this selection will be chosen. Any lineage such as species/kingdom/phylum can be utilized as long as it is contained within the Taxonomic Database. If it is not located within the database, EnTAP will stop the execution immediately and let you know! 

This feature can be utilized with the |flag_taxon| flag. An example command utilizing both common contaminants and a species taxon can be seen below:

.. code-block:: bash

    EnTAP --runN -i path/to/transcriptome.fasta -d path/to/database.dmnd -c fungi -c bacteria --taxon sapiens


.. _over-label:

Picking Up Where You Left Off
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    * annotation_results.emapper.annotations
    * annotation_results_no_hits.emapper.annotations


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
#. Best Hit Selection
#. Gene Ontology / Gene Families

With this functionality of EnTAP, you can execute whatever states you would like with certain commands. Using a '+' will execute from that state to the end, while using a 'x' will stop at that state. These basic commands can be combined to execute whatever you would like. It's easier if I lay out some examples:

* (- - state 1+)
    * This will start at expression filtering and continue to the end of the pipeline

* (- - state 1+5x)
    * This will start at expression filtering and stop at best hit selection

* (- - state 4x)
    * This will just execute similarity search and stop

* (- - state 1+3x5+)
    * This will essentially execute every stage besides similarity searching
    * This is an example of something that may fail, since best hit selection relies on similarity search results

The default 'state' of EnTAP is merely '+'. This executes every stage of the pipeline (or attempts to if the correct commands are in place). 
