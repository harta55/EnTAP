Expression Analysis Flags
=============================

These are the flags specific to Expression Analysis using RSEM. These will be used via the command line (denoted CMD) or ini file (denoted INI).

-a / *-*-align [string] [CMD]
--------------------------------
* Path to alignment file (either SAM or BAM format)
* **Note:** Ignoring this flag will skip expression filtering
* If you have ran alignment with single end reads be sure to use the - -single-end flag as well (paired-end is default)
* Be sure to specify an FPKM threshold

*-*-fpkm [decimal] [INI]
-------------------------------
    * Specify FPKM cutoff for expression filtering
    * Default: 0.5

*-*-single-end [INI]
------------------------
    * Signify your reads are single end for RSEM execution
    * Default: paired-end 