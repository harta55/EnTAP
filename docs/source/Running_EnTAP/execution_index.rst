Executing EnTAP
====================

Execution is the main stage of EnTAP that will annotate a transcriptome input by the user. After following the installation/configuration steps at least once, all of the required databases have been downloaded and configured so that Execution can be ran. Configuration will not need to be ran again unless you would like to update any databases or configure more.

The following stages will be run:

.. toctree::
   :maxdepth: 1
   :numbered:

   Expression Analysis (optional)<Execution/execution_exp_analysis.rst>
   Frame Selection (optional) <Execution/execution_frame_select.rst>
   Similarity Searching <Execution/execution_sim_search.rst>
   Gene Family/Ontology Analysis <Execution/execution_ontology.rst>
   Horizontal Gene Transfer (optional) <Execution/execution_hgt.rst>


Input Files:
^^^^^^^^^^^^^^^^^
Required:

* .FASTA formatted transcriptome file (either protein or nucleotide)
* .dmnd (DIAMOND) indexed databases, which can be formatted in Configuration

Optional:

* .BAM/.SAM alignment file. If left unspecified expression filtering will not be performed. Refer to Expression Analysis section for further detail
* GFF, Donor, and Recipient databases. If left unspecified Horizontal Gene Transfer analysis will not be performed

Sample Run:
^^^^^^^^^^^^^^^^^

A specific run flag (**runP/runN**) must be used:

* runP: Indicates blastp. Frame selection will be ran if nucleotide sequences are input
* runN: Indicates blastx. Frame selection will not be ran with this input


An example run with a nucleotide transcriptome (transcriptome.fnn), two reference DIAMOND databases, an alignment file (alignment.sam), and 8 threads. Expression analysis, frame selection, similarity search, and EggNOG analysis will be ran with the following.

In |run_ini_file_format|, change the following lines for our frame selection test:

.. code-block:: bash

    input=path/to/transcriptome.fnn
    database=path/to/database.dmnd,path/to/database2.dmnd
    align=path/to/alignment.sam

.. code-block:: bash

    EnTAP --runP --run-ini path/to/|config_file| --entap-ini path/to/|run_ini_file|

Further detail on each stage of EnTAP can be found in the page for each stage!

Resuming an EnTAP Run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to save time and make it easier to do different analyses of data, EnTAP allows for picking up where you left off if certain stages were already ran and you'd like to analyze data with different contaminant flags or taxonomic favoring (or more!). As an example, if similarity searching was ran previously you can skip aligning against the database and just analyze the data to save time. This is done by default, if you'd not like to do this set the --overwrite flag to TRUE within the |run_ini_file_format| file.  

In order to pick up and skip re-running certain stages again, the files that were ran previously **must** be in the same directories and have the same names. With an input transcriptome name of 'transcriptome' and example DIAMOND database of 'complete.protein.dmnd':

* Expression Filtering
    * transcriptome.genes.results

* Frame Selection
    * transcriptome.fasta.faa
    * transcriptome.fasta.fnn
    * transcriptome.fasta.lst

* Similarity Search
    * blastp_transcriptome_complete.protein.out

* Gene Family (EggNOG)
    * blastp_transcriptome_eggnog_proteins.out (for runP)
    * blastp_transcriptome_eggnog_proteins.out (for runN)
	
* Horizontal Gene Transfer
   * blastp_transcriptome_complete.protein.out


In order to resume a run, some things must stay the same! Deviations from these flags will cause those stages to be ran again:

* runN/runP

* ontology

* protein

* input

* align

* database

    * Does not necessarily need to remain the same. If additional databases are added, EnTAP will recognize the new ones and run similarity searching on them whilst skipping those that have already been ran

* no-trim

* out-dir

* hgt-donor/hgt_recipient

    * Similar to the DIAMOND databases above, if additional databases are added EnTAP will recognize the new ones and run similarity searching on them whilst skipping those that have already been ran
	
* hgt-gff