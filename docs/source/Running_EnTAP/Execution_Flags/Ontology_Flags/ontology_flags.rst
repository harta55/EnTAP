Ontology Flags
=====================

These are flags specific to the Ontology portion of the pipeline using either EggNOG or InterProScan. These will be used via the command line (denoted CMD) or ini file (denoted INI).

General Flags
------------------

*-*-ontology [multi-integer] [INI]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    * Specify which ontology packages you would like to use

        * 0 - EggNOG (default)
        * 1 - InterProScan

    * Both or either can be specified with multiple flags

        * Ex: - - ontology 0 - - ontology 1
        * This will run both EggNOG and InterProScan 

InterProScan Specific Flags
------------------------------------------

*-*-protein [multi-string] [INI]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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