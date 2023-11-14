Changelog
==================
This page contains (mostly) all of the changes that were made between each version of EnTAP. 

EnTAP v1.0.1 (November 13, 2023)
------------------------------------------
    * Fixed an issue with the formatting of the '...gene_ontology_terms.tsv' file output from EnTAP 

EnTAP v1.0.0 (September 26, 2023)
------------------------------------------
    * Updated RSEM (v1.3.3), TransDecoder (v5.7.1), and DIAMOND (v2.1.8) libraries in EnTAP repository. EnTAP is now compatible with these versions
    * Added new test data (under 'test_data') directory in the EnTAP repository. This should work with latest versions of software being used by EnTAP. See docs for how to run.
    * Added additional statistics/percentages at the end of the Log File
    * Added support for Tidyverse format (TSV's will now print 'NA' for empty data)
    * Added Dockerfile to EnTAP repository
    * Added support for a new Gene Ontology term TSV output format. Similar to other formats, but combined into one file. More info can be seen with '--output--format' flag
    * Changed DIAMOND runs (during Ontology and during Similarity Searching) to use 'very-sensitive' from 'more-sensitive'. This should give more alignments, but may take longer to execute now
    * Removed the '--level' flag for Gene Ontology levels. This was not useful to users and caused confusion. Instead, all GO Terms will be printed by default and gene ontology levels are removed from output
    * Renamed and restructured many of the output in the 'final_results' directory for better clarity
    * Fixed issue if trying to run EnTAP configuration locally to rebuild the EnTAP database. Users may have seen it fail during the Gene Ontology stage due to a change in formatting of the Gene Ontology database

EnTAP Beta v0.10.9-beta (June 29, 2023)
------------------------------------------
    * Optimizations made to similarity searching using the EnTAP SQL database. Should improve speed
    * Added 'api-taxon' command that will verify whether an input taxon can be found in the taxonomy database. It will return json formatted text
    * Added additional messaging throughout EnTAP execution to stdout
    * Removed Test Data from repository, it is no longer compatible with latest version of software within EnTAP. Will add back an updated dataset next version
    * Changed DIAMOND command to use '--max-target-seqs' instead of '--top' command
    * Fixed an issue where duplicate sequences were printed to the final_annotations files
    * Fixed an issue where the taxanomic species may not have been found when searching against the SQL EnTAP database

EnTAP Beta v0.10.8-beta (March 21, 2021)
------------------------------------------
    * This version requires a new version of the EnTAP database to be downloaded
    * Added Gene Enrichment files as an output option(gene ID + effective length and geneID + GO term). These can be seen with the output-type flag in the ini file
    * Changed Gene Ontology level printing. 0 will continue to print every term. Other levels will now print that level AND higher. So a level of 1 will print 1, 2, 3, etc. Previous a level of 1 would only print GO Terms with a level of 1
    * Changed 'uninformative' input from a file to a list of terms in the ini file. Much more straightforward this way
    * If no alignments are found against a database during DIAMOND, the pipeline will no longer exit, it will continue to the next database. If no alignments are found against any databases, it will stop at that point
    * Fixed a bug where TransDecoder output may not have been parsed correctly for some users. This presented itself as a parsing error and halted EnTAP at that stage of the pipeline
    * Fixed bug where InterProScan Mobidlite database was giving an error for some users (and halting execution)

EnTAP Beta v0.10.7-beta (October 6, 2020)
------------------------------------------

    * Fixed an issue where certain sequence headers may not have been parsed properly resulting in unrecognized sequence errors during Similarity Searching

EnTAP Beta v0.10.6-beta (August 26, 2020)
------------------------------------------

    * Added support to pipe the TransDecoder flag '--no_refine_starts' during Execution
    * Fixed an issue where error messages during EggNOG searching would not get printed (seg fault)
    * Contaminant information will not be printed to the log if there are none

EnTAP Beta v0.10.5-beta (August 12, 2020)
------------------------------------------

    * Added a step to remove the stop codon ('*') sometimes printed at the end of the TransDecoder FASTA output. This may have caused an issue when running TransDecoder and InterProScan together

EnTAP Beta v0.10.4-beta (July 29, 2020)
------------------------------------------

    * Fixed an issue where expression analysis transcriptome generation would fail (error message presented to user as 'frame selection')

EnTAP Beta v0.10.3-beta (July 28, 2020)
------------------------------------------

    * Fixed a parsing issue of user inputs for contanminants and taxon

EnTAP Beta v0.10.2-beta (July 26, 2020)
------------------------------------------

    * Fixed a pathing issue when EnTAP generated frame selected transcriptomes

EnTAP Beta v0.10.1-beta (July 19, 2020)
------------------------------------------

Note: Please use v0.10.2-beta or later instead of this version

    * Added support for TransDecoder for Frame Selection
    * Added TPM as an additional output from Expression Filtering
    * Added an .ini file and moved many commands/paths from the command line to this
    * Standardized/finalized output header namings for gFACs support
    * Changed the default Frame Selection software to TransDecoder. GeneMarkS-T can still be selected through the .ini file
    * Changed the default Gene Ontology level to 1. This can be easily changed through the ini file
    * Fixed issue where some EggNOG descriptions were not printed to the final output
    * Fixed a few issues with older GCC versions
    * Fixed an issue where GeneMarkS-T would write to the working directory

EnTAP Beta v0.9.2-beta (June 4, 2020)
------------------------------------------

    * Updated EggNOG Database links


EnTAP Beta v0.9.1-beta (January 12, 2020)
-------------------------------------------

    * Changed --trim flag to --no-trim. Trimming sequence headers to the first space is the default now. If you have executions from previous versions, you may need to use the --no-trim flag as needed for backwards compatibility (picking up where you left off)
    * Fixed a bug where the --single-end command was not properly recognized

EnTAP Beta v0.9.0-beta (May 12, 2019)
-------------------------------------------

    * This release focused on reducing installation complexity and removing dependencies
    * Overhauled the configuration/execution process by removing EggNOG-mapper and replacing it with an internal EnTAP method. This will make installation and both stages much clearer for the user
    * Removed Boost Libraries from dependencies further reducing installation complexity
    * Added printing of error messages to the standard log from any software being used by EnTAP. This will make debugging much easier
    * Added UniProt mapping to the EnTAP database. This will pull any additional mapping information from UniProt Swiss-Prot alignments
    * Updated supported DIAMOND version to 0.9.9
    * The EnTAP database MUST be re-configured for this release
    * Resolved any incompatibility with DIAMOND and EggNOG databases as well as versioning problems
    * Standardized EnTAP log entries and added additional statistics
    * - -ontology flag will now use EnTAP's method of EggNOG accession (0) or InterProScan (1)
    * Bug fixes


EnTAP Beta v0.8.4-beta (August 2, 2018)
------------------------------------------------

    * Fixed an issue when inputting already translated sequences


EnTAP Beta v0.8.3-beta (May 23, 2018)
------------------------------------------

    * Minor bug fixes
    * Changes to CMake to hopefully resolve issues a couple users had with linking to Boost Libraries


EnTAP Beta v0.8.2-beta (April 29, 2018)
-------------------------------------------

    * Revamped configuration stage of EnTAP (reduced time and hopefully made things clear/more compatible across systems)
    * Removed - -database-out flag (seemed a bit redundant to me). - -outfiles flag will be the default when indexing databases
    * Added - -data-generate flag. This can be specified in EnTAP config stage (no effect during execution) for whether you'd like to generate the EnTAP databases rather than downloading from FTP address
    * Added - -data-type flag. This can be used in either configuration or execution. Specifies which database you'd like to download/generate or use during execution. Binary (0, default) or SQL (1). Binary is faster with more memory usage, SQL will be slower but easier compatibility.
    * Combined EnTAP databases into one (entap_database.sql/entap_database.bin). WARNING: Re-download or configuration of databases is REQUIRED with this newer version.
    * Removed download_tax.py script (no longer necessary)


EnTAP Beta v0.8.1-beta (April 14, 2018)
------------------------------------------

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


EnTAP Beta v0.8.0-beta (December 16, 2017)
-------------------------------------------------

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
		

EnTAP Beta v0.7.4.1-beta (September 5, 2017)
--------------------------------------------------

    * Minor changes to taxonomic database download and indexing

EnTAP Beta v0.7.4-beta (August 26, 2017)
----------------------------------------------

    * Initial beta release!
