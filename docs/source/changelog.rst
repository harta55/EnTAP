Changelog
==================
This page contains (mostly) all of the changes that were made between each version of EnTAP. The current latest version is EnTAP Beta v0.8.4-beta

EnTAP Beta v0.8.4-beta
-----------------------
    * Fixed some protein input issues

EnTAP Beta v0.8.3-beta
------------------------

    * Minor bug fixes
    * Changes to CMake to hopefully resolve issues a couple users had with linking to Boost Libraries


EnTAP Beta v0.8.2-beta
------------------------

    * Revamped configuration stage of EnTAP (reduced time and hopefully made things clear/more compatible across systems)
    * Removed - -database-out flag (seemed a bit redundant to me). - -outfiles flag will be the default when indexing databases
    * Added - -data-generate flag. This can be specified in EnTAP config stage (no effect during execution) for whether you'd like to generate the EnTAP databases rather than downloading from FTP address
    * Added - -data-type flag. This can be used in either configuration or execution. Specifies which database you'd like to download/generate or use during execution. Binary (0, default) or SQL (1). Binary is faster with more memory usage, SQL will be slower but easier compatibility.
    * Combined EnTAP databases into one (entap_database.sql/entap_database.bin). WARNING: Re-download or configuration of databases is REQUIRED with this newer version.
    * Removed download_tax.py script (no longer necessary)


EnTAP Beta v0.8.1-beta
------------------------

    * Added additional error logging to provide more information when something goes wrong
    * Configuration file mandatory (default place to look is current working directory)
    * Changed tax database paths in config file to avoid confusion (separate text and bin). Config file must be re-downloaded/generated!
    * Defaults/output during configuration changed to config file then if not found, database-out flag
    * Added deletion of empty files if a certain stage failed (preventing re-reading an empty file)
    * Added errors/warnings for no alignments/hits in each stage
    * entap_out directory changed to transcriptomes to be more clear (holds only transcriptomic data)
    * Final EnTAP output files moved from the root outfiles directory to final_results directory
    * Several filename changes to add consistency in new transcriptomes directory (final transcriptome is now _final.fasta. 
    * Several title changes to the log file to mitigate confusion
    * EggNOG no longer broken down into separate files - those that hit and those that did not hit a database. Now entire transcriptome is pushed with one output file
    * 10 species/contaminants/other in similarity searching statistics has been changed to 20 to provide more information to the user
    * Best hit selection state combined with similarity search
    * Added 'N' as an accepted nucleotide
    * Several behind the scenes changes
    * Fixed Cmake global installation issue
    * Fixed incorrect error codes
    * Fixed InterPro printing bug to no hits/hits files
    * Fixed Frame Selection not printing new lines for certain files


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
