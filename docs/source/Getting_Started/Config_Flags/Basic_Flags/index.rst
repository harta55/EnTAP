Basic Flags
=====================

These are the general flags that can be used during Configuration.

-d / *-*-database [string]
-------------------------------
* Specify any number of FASTA formatted databases you would like to configure for EnTAP
* Not necessary if you already have DIAMOND configured databases (.dmnd)

*-*-out-dir [string]
----------------
* Specify an output directory for the databases to be sent to (recommended)
* This will send the EnTAP database and DIAMOND databases to this location

-t / *-*-threads [number]
---------------------
* Specify thread number for Configuration

*-*-data-generate
------------------------
* Specify this flag is you would like to generate the EnTAP database rather than downloading from FTP (default)
* I'd only use this if you're having issues with the FTP

*-*-data-type [flag]
-------------------
* Specify which databases you'd like to generate/download

    * 0. Binary Database (default) - This will be much quicker and is recommended
    * 1. SQL Database - Slower although will be more easily compatible with every system

* This can be flagged multiple times (ex: - - data-type 0 - - data-type 1)
* I would not use this flag unless you are experiencing issues with the EnTAP Binary Database