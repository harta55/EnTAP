Similarity Search Flags
=============================

These are the flags specific to Similarity Searching using DIAMOND.

-d / *-*-database
-----------------------
* Specify any number of FASTA formatted databases you would like to configure for EnTAP
* Not necessary if you already have DIAMOND configured databases (.dmnd)

*-*-data-type
-------------------
* Specify which EnTAP database you'd like to use for execution (UniProt, Gene Ontology, and Taxonomy lookups)

    * 0. Binary Database (default) - This will be much quicker and is recommended
    * 1. SQL Database - Slower although will be more easily compatible with every system

* This can be flagged multiple times (ex: - - data-type 0 - - data-type 1)
* I would not use this flag unless you are experiencing issues with the EnTAP Binary Database

-c / *-*-contam [multi-string]
----------------------------------------
    * Specify :ref:`contaminant<tax-label>` level of filtering
    * Multiple contaminants can be selected through repeated flags

*-*-taxon [string]
------------------
    * This flag will allow for :ref:`taxonomic<tax-label>` 'favoring' of hits that are closer to your target species or lineage. Any lineage can be used as referenced by the NCBI Taxonomic database, such as genus, phylum, or species.
    * Format **must** replace all spaces with underscores ('_') as follows: "- -taxon homo_sapiens" or "- -taxon primates"

*-*-level [multi-integer]
--------------------------------
    * Specify Gene Ontology levels you would like to normalize to
    * Any amount of these flags can be used
    * Default: 0 (every level), 3, 4
    * More information at: http://geneontology.org/page/ontology-structure

-e [scientific]
-----------------
    * Specify minimum E-value cutoff for similarity searching (scientific notation)
    * Default: 10E-5

*-*-tcoverage [decimal]
-----------------------------
    * Specify minimum target coverage for similarity searching
    * Default: 50%

*-*-qcoverage [decimal]
--------------------------
    * Specify minimum query coverage for similarity searching
    * Default: 50%

*-*-uninformative [string]
----------------------------------
    * Path to a list of terms you would like to be deemed "uninformative"
    * The file **must** be formatted with one term on each line of the file
    * Example (defaults):
    
        * conserved
        * predicted
        * unnamed
        * hypothetical
        * putative
        * unidentified
        * uncharacterized
        * unknown
        * uncultured
        * uninformative