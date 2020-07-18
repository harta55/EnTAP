Required Flags
=====================

These are the required flags that must be used during Execution.

*-*-runP / *-*-runN
-------------------------
* Specify a blastp or blastx annotation
* If - -runP is selected with a nucleotide input, frame selection will be ran and annotation stages will be executed with protein sequences (blastp)
* If - -runP is selected with a protein input, frame selection will not be ran and annotation will be executed with protein sequences (blastp)
* If - -runN is selected with nucleotide input, frame selection will not be ran and annotation will be executed with nucleotide sequences (blastx)

-i / *-*-input [string]
-------------------------
* Path to the transcriptome file (either nucleotide or protein)

-d / *-*-database [multi-string]
-------------------------
    * Specify up to 5 DIAMOND indexed (.dmnd) databases to run similarity search against

*-*-|flag_path| [string]
---------------------
* Point to |config_file| to specify file paths
* Default: |config_file| residing in the current working directory