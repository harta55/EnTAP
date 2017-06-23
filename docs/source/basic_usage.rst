Basic Usage
============

*enTAP* has two stages of execution, :ref:`configuration<config-label>` and :ref:`run<run-label>`. Configuration is generally ran first (and may only need to be ran once) to setup databases while run is reserved for the main annotation pipeline.

.. _config-label:

Configuration
-------------
Configuration is the first stage of *enTAP* that will download and configure the necessary databases for full functionality. Additionally, databases in traditional .fasta (or similar) format must be configured to run with *enTAP* for faster similarity searching. This can be done with any database you have previously downloaded and will be configured and sent to the /bin folder within the *enTAP* directory. 

To run configuration with a sample database, the command is as follows:

.. code-block:: bash

    enTAP --config -d path/to/database

This stage must be done at least once prior to :ref:`running<run-label>`. Once the database is configured, you need not do it again unless you updated your original database or plan on configuring several others.


**Note:** If you already have DIAMOND (.dmnd) configured databases, you can skip the configuration of that database. Although, due to other database downloading (taxonomy and ontology), you must still run the configuration stage without any flags


**Note:** This is the only stage that requires connection to the internet.

Flags:
^^^^^^^^^^^^^^^^^^^^^

Required Flags:

* The only required flag is **--config**. Although in order to run the full *enTAP* pipeline, you must have a .dmnd configured database.


Optional Flags:

* -d : Specify any number of databases you would like to configure for *enTAP*


Memory Usage:
^^^^^^^^^^^^^^

Memory usage will vary depending on the number of database you would like configured. Although, *enTAP* will download several other databases as well:

* Gene Ontology References: 6Mb
* NCBI Taxonomy: 400Mb

.. _run-label:

Run
-------------
The run stage of *enTAP* is the main annotation pipeline. After configuration is ran at least once, this can be ran continually without requiring configuration to be ran again. 

Input Files:
^^^^^^^^^^^^
Required:

* .FASTA formatted transcriptome file (either protein or nucleotide)
* .dmnd (DIAMOND) indexed databases, which can be formatted in the :ref:`configuration<config-label>` stage. Up to 4 can be chosen


Optional:

* .BAM/.SAM alignment file. If left unspecified expression filtering will not be performed. 

Sample Run:
^^^^^^^^^^^

A specific run flag (**runP/runN**) must be used:

* runP: Indicates protein input transcripts. **Note:** Selection of this option will skip the frame selection portion of the pipeline.
* runN: Indicates nucleotide input transcripts. Selection of this option will cause frame selection to be ran. 


An example run with a nucleotide transcriptome:

.. code-block:: bash

    enTAP --runN -i path/to/transcriptome.fasta -d path/to/database.dmnd -d path/to/database2.dmnd -a path/to/alignment.sam


With the above command, the entire *enTAP* pipeline will be ran. However, it is entirely possible to skip frame selection by inputting protein transcripts (--runP) or skip expression filtering by excluding an alignment file. 


Flags:
^^^^^^^^^^^^^^^^^^^^^

Required Flags:

* (--runP/--runN)
    * Specification of input transcriptome file. runP for protein (skip frame selection) or runN for nucleotide (frame selection will be ran)

* (-i/--input)
    * Path to the transcriptome file (either nucleotide or protein)

* (-d/--database)
    * Specify up to 4 DIAMOND indexed (.dmnd) databases to run similarity search against

Optional Flags:

* (-a/--align)
    * Path to alignment file (either SAM or BAM format)
    * **Note:** Ignoring this flag will skip expression filtering

* (--contam)
    * Specify :ref:`contaminant<tax-label>` level of filtering
    * Multiple contaminants can be selected through repeated flags

* (--species)
    * Species of your transcriptome. This will allow for taxonomic 'favoring' of hits that are closer to your target species.
    * Format **must** be as follows: "--species homo_sapiens"

* (--tag)
    * Specify output folder labelling.
    * (default: outfiles)


.. _tax-label:

Taxonomic Contaminant Filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


