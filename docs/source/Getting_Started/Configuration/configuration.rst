Configuration
=====================

After Installation is complete, EnTAP must be configured for use. This stage will simply download and configure the necessary databases for full functionality and only needs to be ran once (unless more DIAMOND databases want to be updated or configured for similarity searching). It is recommended to run this once after installation, then users of EnTAP can run it additional times if more databases are required. Before starting Configuration, make sure to review both EnTAP ini files and ensure at least the DIAMOND and database flags are setup properly. 

This will be broken up into two sections. If you would like to configure databases for searching with DIAMOND, refer to :ref:`Preparing Your Reference DIAMOND Databases<usage-label>`. If you already have DIAMOND databases configured, simply move on to :ref:`Running Configuration<usage_config-label>`. 

.. _usage-label:

Preparing Your Reference DIAMOND Databases
----------------------------------------------

All source databases must be provided in FASTA format (protein) so that they can be indexed for use by DIAMOND.  This can be completed independent of EnTAP with DIAMOND (- - makedb flag) or as part of the Configuration phase of EnTAP. This section will focus on downloading and preparing some of the more common FASTA source databases. If you already have DIAMOND databases configured, you can skip to :ref:`Running Configuration<usage_config-label>`. Even if you have a DIAMOND database already configured, Configuration must still be ran!

While any protein FASTA database can be used, it is recommended to use NCBI (Genbank) sourced databases such as RefSeq databases or NR.  In addition, EnTAP can easily accept EBI databases such as UniProt/SwissProt.  

EnTAP can recognize the species information from these header formats ONLY (NCBI and UniProt):

* [homo sapiens]

* OS=homo sapiens

If the individual FASTAs in a custom database you create do not adhere to one of these two formats, it will not be possible to weight taxonomic or contaminant status from them. You will need to change the headers to ensure they align. 

The following FTP sites contain common reference databases that EnTAP can recognize:
   * RefSeq: |ref_comp|

   * Plant RefSeq: |ref_plant|

   * Mammalian Vertebrate RefSeq: |ref_mamm|

   * Other Vertebrate RefSeq: |ref_vert_other|

   * Invertebrate RefSeq: |ref_invert|

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

Running Configuration
-------------------------------

As mentioned before, Configuration will only need to be ran once to download necessary files for EnTAP unless you would like to configure/update more databases. 

In general, the command to run EnTAP Configuration is as follows. Before execution, please refer to the relevant Configuration Flags before running this to be sure it goes smoothly and duplicate databases are not downloaded:

.. code-block:: bash

    EnTAP --config --run-ini path/to/entap_config.ini --entap-ini path/to/entap_run.params

To run configuration and format a FASTA database for DIAMOND, simply add these databases to the |run_ini_file_format| file.

.. warning ::
    Sometimes DIAMOND database versions are not always cross-compatible with different versions of DIAMOND. To avoid this, configure databases and eventually execute with the same version of DIAMOND

In both cases, the following databases will be downloaded and configured:

* EnTAP Database:
    * Comprised of Gene Ontology, UniProt, and Taxonomic mappings for use during Execution. FTP downloaded file.
    * Downloaded from |entap_bin_ftp|
    * Filename: entap_database.bin
    * The SQL version is the same database, but formatted as a SQL database. Only one version of the database is needed (binary is used by default and SQL is much slower but uses less memory)
    * If you experience any trouble in downloading, you can simply specify the - - data-generate flag during Configuration to configure it locally (more on that later)
    * The database for the newest version of EnTAP will always reside in the "latest" FTP directory. Keep in mind, if you are using an older version of EnTAP, you do not want to download from the "latest" directory. Instead, you will need to consider the version you are using. The FTP will always be updated only when a new database version is created. For example, if you see v0.8.2 and v0.8.5 on the FTP while you are using v0.8.3, you will download the database located in the v0.8.2 directory. 

* EggNOG DIAMOND Reference:
    * Reference database containing EggNOG database entries
    * DIAMOND formatted database is downloaded from |eggnog_dmnd_ftp|
    * Filename: eggnog_proteins.dmnd

* EggNOG SQL Database:
    * SQL database containing EggNOG mappings
    * Downloaded from |eggnog_sql_ftp|
    * Filename: eggnog.db
    * Note, when referencing this file in the |config_file_format|, you must use the directory that contains this file with the --egg-map-data flag, rather than the path to the file itself

.. note:: Either the EnTAP binary database (default) or the EnTAP SQL database is required for execution. Both are not needed.