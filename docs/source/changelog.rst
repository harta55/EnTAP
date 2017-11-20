Changelog
==================
This page contains (mostly) all of the changes that were made between each version of EnTAP. The current latest version is EnTAP Beta v0.8.0-beta

EnTAP Beta v0.8.0-beta
------------------------

    * Overhaul of the taxonomic/gene ontology databases
        
        * Faster accession/indexing
        * MUST be re-downloaded and re-indexed (or use the updated versions that come with the EnTAP distribution)
        * Taxonomic database includes thousands more entries with synonyms of many species
        * Perl is no longer a dependency, with Python being used to download the database

    * Added blastx support

        * Blastx now allowed for ALL stages of annotation (similarity search + ontology)
        * --runN flag now specifies blastx (frame selection will not be ran)
        * --runP flag now specifies blastp (frame selection will be performed if nucleotide sequences are input)
        
    * Added InterProScan support

        * Now possible to run EggNOG and/or InterProScan (with both blastx or blastp)
        * EggNOG and/or InterProScan specified with --ontology flag (0 and/or 1)
        * Full output of both will be provided in the final annotations file
        
    * Added additional statistics to the log file for EggNOG and Expression Analysis
    * Added numerous file/path/software checks to the start of an EnTAP run

        * Test runs/path checks are performed on all software that will be ran
        * Additional checks to specific flags
        * These checks can be turned off for an EnTAP run with --no-check flag (not advised!) 

    * --tag flag changed to --out-dir to specify output directory (not just what you'd like it named as)
  
        * Defaults to current directory with /outfiles folder

    * --paired-end flag for Expression Filtering changed to --single-end (with paired-end being the default)
    * Added contaminant and informative yes/no columns in final annotations file (among other headers)
    * Added ability to input your own list of informative/uninformative terms for EnTAP to flag
    * Added contaminant and none contaminant final annotation files
    * Fixed a sequence id issue in Expression Filtering not mapping to BAM/SAM file
    * Fixed a bug in --trim flag for sequence headers
    * Fixed a bug where some systems had issues with graphing
    * Debug and log files are now time stamped and not overwritten
    * Fixed pathing for EnTAP configuration and made more streamlined
    * Fixed several instances of older compilers complaining
    * Added a lot of error messaging to help diagnose any issues easily
    * Changed similarity search to have full database name, not path
    * Fixed a bug in parsing input fasta file (added corrupt file checks)
		

EnTAP Beta v0.7.4.1-beta
-------------------

    * Minor changes to taxonomic database download and indexing

EnTAP Beta v0.7.4-beta
------------------

    * Initial beta release!
