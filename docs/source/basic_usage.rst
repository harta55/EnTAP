.. _NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/taxonomy
.. |libs_dir| replace:: /libs
.. |entap_dir| replace:: /EnTAP
.. |src_dir| replace:: /src
.. |config_file| replace:: entap_config.txt
.. |bin_dir| replace:: /bin
.. |data_dir| replace:: /databases
.. |tax_file| replace:: download_tax.pl
.. |graph_file| replace:: entap_graphing.py
.. |go_term| replace:: go_term.entp
.. |tax_bin| replace:: ncbi_tax_bin.entp
.. |tax_data| replace:: ncbi_tax.entp

.. |flag_path| replace:: paths


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
    * genemarkst_exe_path=/EnTAP//libs/gmst_linux_64/gmst.pl
    * eggnog_exe_path=/EnTAP/libs/eggnog-mapper/emapper.py
    * eggnog_download_exe=/EnTAP/libs/eggnog-mapper/download_eggnog_data.py
    * eggnog_database=/EnTAP/libs/eggnog-mapper/data/eggnog.db (downloaded during Configuration)
    * entap_tax_database=/EnTAP/bin/ncbi_tax_bin.entp (binary version, downloaded during Configuration)
    * entap_tax_download_script=/EnTAP/src/download_tax.pl
    * entap_go_database=/EnTAP/bin/go_term.entp (binary version, downloaded during Configuration)
    * entap_graphing_script=/EnTAP/src/entap_graphing.py


These can be changed to whichever path you would prefer. If something is globally installed, just put a space " " after the '='. EnTAP will recognize these paths first and they will override defaults. 


This configuration file will be automatically detected if it is in the same directory as the EnTAP .exe, otherwise the path to it can be specified through the |flag_path| flag. 

.. note:: Be sure you set the paths before moving on (besides the databases that haven't been downloaded yet)!


.. _usage-label:

Usage
^^^^^

All source databases must be provided in FASTA format so that they can be indexed for use by DIAMOND.  This can be completed independent of EnTAP with DIAMOND or as part of the configuration phase of EnTAP.  While any FASTA database can be used, it is recommended to use NCBI (Genbank) sourced databases such as RefSeq databases or NR.  In addition, EnTAP can easily accept EBI databases such as UniProt/SwissProt.  EnTAP can read the species information from these header formats.  If the individual FASTAs in a custom database do not adhere to one of these two formats, it will just not be possible to weight examine taxanomic or contaminant status from them.  

The following FTP sites contain common reference databases that enTAP can recognize:
   * RefSeq:
   * Arthropod RefSeq:
   * Plant RefSeq:
   * Mammalian RefSeq:
   * NR:
   * SwissProt:
   * UniProt:
....

It is generally recommended that a user select at least three databases with varying levels of NCBI curation.  Unless the species is very non-model (i.e. does not have close relatives in databases such as RefSeq, it is not necessary to use the full NR database which is less curated).


To run configuration with a sample database, the command is as follows:

.. code-block:: bash

    EnTAP --config -d path/to/database

This stage must be done at least once prior to :ref:`running<run-label>`. Once the database is configured, you need not do it again unless you updated your original database or plan on configuring several others.


.. note:: If you already have DIAMOND (.dmnd) configured databases, you can skip the configuration of that database. Although, due to other EnTAP database downloading (taxonomy and ontology), configuration must still be ran at least once without any flags.

Configuration can be ran without formatting a database as follows:

.. code-block:: bash

    EnTAP --config


.. note:: This is the only stage that requires connection to the Internet.

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


Memory Usage:
^^^^^^^^^^^^^^

Memory usage will vary depending on the number of databases you would like configured. Although, EnTAP will download several other databases as well:

* Gene Ontology References: 6Mb
* NCBI Taxonomy: 400Mb
* EggNOG Database: 30Gb

....

.. _run-label:

Run
-------------
The run stage of *EnTAP* is the main annotation pipeline. After configuration is ran at least once, this can be ran continually without requiring configuration to be ran again (unless more databases will be configured). 

Input Files:
^^^^^^^^^^^^
Required:

* .FASTA formatted transcriptome file (either protein or nucleotide)
* .dmnd (DIAMOND) indexed databases, which can be formatted in the :ref:`configuration<config-label>`stage. 

Optional:

* .BAM/.SAM alignment file. If left unspecified expression filtering will not be performed. 

Sample Run:
^^^^^^^^^^^

A specific run flag (**runP/runN**) must be used:

* runP: Indicates protein input transcripts. Selection of this option will skip the frame selection portion of the pipeline.
* runN: Indicates nucleotide input transcripts. Selection of this option will cause frame selection to be ran. 


An example run with a nucleotide transcriptome:

.. code-block:: bash

    EnTAP --runN -i path/to/transcriptome.fasta -d path/to/database.dmnd -d path/to/database2.dmnd -a path/to/alignment.sam


With the above command, the entire EnTAP pipeline will run. Both frame selection and expression filtering can be skipped if preferred by the user.  EnTAP would require protein sequences (indicated by --runP) in order to avoid frame selection.  If there is not a short read alignment file provided in SAM/BAM format, then expression filtering via RSEM will be skipped. 


Flags:
^^^^^^^^^^^^^^^^^^^^^

Required Flags:

* (- -runP/- -runN)
    * Specification of input transcriptome file. runP for protein (skip frame selection) or runN for nucleotide (frame selection will be ran)

* (-i/- -input)
    * Path to the transcriptome file (either nucleotide or protein)

* (-d/- -database)
    * Specify up to 4 DIAMOND indexed (.dmnd) databases to run similarity search against

Optional Flags:

* (-a/- -align)
    * Path to alignment file (either SAM or BAM format)
    * **Note:** Ignoring this flag will skip expression filtering
    * If you have ran alignment with paired end reads be sure to use the - -paired-end flag as well

* (- - contam)
    * Specify :ref:`contaminant<tax-label>` level of filtering
    * Multiple contaminants can be selected through repeated flags

* (- - species)
    * This flag will allow for :ref:`taxonomic<tax-label>` 'favoring' of hits that are closer to your target species or lineage. Any lineage can be used as referenced by the NCBI Taxonomic database, such as genus, phylum, or species.
    * Format **must** replace all spaces with underscores ('_') as follows: "- -species homo_sapiens" or "- -species primates"

* (- - level)
    * Specify Gene Ontology levels you would like to normalize to
    * Any amount of these flags can be used

* (- - tag)
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

* (- - paired-end)
    * Signify your reads are paired end for RSEM execution

* (- - graph)
    * This will check whether or not your system has graphing functionality supported
    * If Python with the Matplotlib module are installed on your system graphing should be enabled!
    * This can be specified on its own

* (-t/ - - threads)
    * Specify the number of threads of execution

* (- - state)
    * Precise control over execution :ref:`stages<state-label>`. This flag allows for certain parts to be ran while skipping others. 
    * Warning: This may cause issues depending on what you plan on running! 


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

During best hit selection of similarity searched results, taxonomic consideration can utilized. If a certain lineage (such as sapiens) is specified, hits closer in taxonomic lineage to this selection will be chosen. Any lineage such as species/kingdom/phylum can be utilized as long as it is contained within the Taxonomic Database


.. _over-label:

Picking Up Where You Left Off
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to save time and make it easier to do different analyses of data, EnTAP allows for picking up where you left off if certain stages were already ran and you'd like analyze data with different contaminant flags or taxonomic favoring. As an example, if similarity searching was ran previously you can skip hitting against the database and analyze the data to save time. However, the - - overwrite flag will not allow for this as it will remove previous runs and not recognize them. 

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
* (-i / - - input)

* (-a / - - align)

* (-d / - - database)
    * Do not necessarily need to remain the same. If additional databases are added, EnTAP will recognize the new ones and run similarity searching on them    

* (- - qcoverage)

* (- - tcoverage)


.. _state-label:

State Control
^^^^^^^^^^^^^^
