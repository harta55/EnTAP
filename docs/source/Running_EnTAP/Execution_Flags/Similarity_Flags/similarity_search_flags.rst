Similarity Search Flags
=============================

These are the flags specific to Similarity Searching using DIAMOND. These will be used via the command line (denoted CMD) or ini file (denoted INI).

-d / *-*-database [CMD]
---------------------------
* Specify any number of FASTA formatted databases you would like to configure for EnTAP
* Not necessary if you already have DIAMOND configured databases (.dmnd)

*-*-data-type [INI]
------------------------
* Specify which EnTAP database you'd like to use for execution (UniProt, Gene Ontology, and Taxonomy lookups)

    * 0. Binary Database (default) - This will be much quicker and is recommended
    * 1. SQL Database - Slower although will be more easily compatible with every system

* This can be flagged multiple times (ex: - - data-type 0 - - data-type 1)
* I would not use this flag unless you are experiencing issues with the EnTAP Binary Database

-c / *-*-contam [multi-string] [INI]
----------------------------------------
    * Specify :ref:`contaminant<tax-label>` level of filtering
    * Multiple contaminants can be selected through repeated flags

*-*-taxon [string] [INI]
-----------------------------
    * This flag will allow for :ref:`taxonomic<tax-label>` 'favoring' of hits that are closer to your target species or lineage. Any lineage can be used as referenced by the NCBI Taxonomic database, such as genus, phylum, or species.
    * Format **must** replace all spaces with underscores ('_') as follows: "- -taxon homo_sapiens" or "- -taxon primates"

*-*-level [multi-integer] [INI]
--------------------------------
    * Specify Gene Ontology levels you would like to normalize to (ex: 0, 1, 2, 3, 4)
    * This should only be used as a rough idea, some of the levels can vary slightly
    * A level of '0' indicates all levels will be printed while a level of 2 will indicate that all levels of 2 or higher will be printed.
    * Any amount of these flags can be used
    * Default: 0, 1
    * More information at: http://geneontology.org/page/ontology-structure

-e [scientific] [INI]
----------------------------
    * Specify minimum E-value cutoff for similarity searching (scientific notation)
    * Default: 10E-5

*-*-tcoverage [decimal] [INI]
----------------------------------
    * Specify minimum target coverage for similarity searching
    * Default: 50%

*-*-qcoverage [decimal] [INI]
----------------------------------
    * Specify minimum query coverage for similarity searching
    * Default: 50%

*-*-uninformative [string] [INI]
----------------------------------
    * Comma-deliminated list of terms you would like to be deemed "uninformative". Any alignments during Similarity Searching tagged as uninformative will be ranked lower.
    * Default: conserved, predicted, unnamed, hypothetical, putative, unidentified, uncharacterized, unknown, uncultured, uninformative