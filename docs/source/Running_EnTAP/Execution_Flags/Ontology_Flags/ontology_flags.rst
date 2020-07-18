Ontology Flags
=====================

These are flags specific to the Ontolgoy portion of the pipeline using either EggNOG or InterProScan.

General Flags
------------------

*-*-level [multi-integer]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    * Specify Gene Ontology levels you would like to normalize to
    * Any amount of these flags can be used
    * Default: 0 (every level), 3, 4
    * More information at: http://geneontology.org/page/ontology-structure

*-*-ontology [multi-integer]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    * Specify which ontology packages you would like to use

        * 0 - EggNOG (default)
        * 1 - InterProScan

    * Both or either can be specified with multiple flags

        * Ex: - - ontology 0 - - ontology 1
        * This will run both EggNOG and InterProScan 

InterProScan Specific Flags
------------------------------------------

*-*-protein [multi-string]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    * Use this option if you would like to run InterProScan
    * Specify databases to run against (you must have them already installed)
      
        * tigrfam
        * sfld
        * prodom
        * hamap
        * pfam
        * smart
        * cdd
        * prositeprofiles
        * prositepatterns
        * superfamily
        * prints
        * panther
        * gene3d
        * pirsf
        * coils
        * mobidblite