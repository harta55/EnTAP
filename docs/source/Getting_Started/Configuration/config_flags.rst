Configuration Flags
============================

These are the flags for the configuration process of EnTAP. These will be used via the command line (denoted CMD) or ini file (denoted INI).

.. list-table:: **Required Flags**
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1
   
   * - param
     - description
     - location (cmd/ini)
     - qualifier
     - example
   * - config
     - Must be specified in order for EnTAP to run Configuration as opposed to normal execution
     - cmd
     - flag
     - config
   * - |flag_path|
     - Point to the |config_file| to specify paths. DIAMOND is only required for Configuration.
     - cmd
     - string
     - /path/to/|config_file|
     
.. list-table:: **Optional Flags**
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1
   
   * - param
     - description
     - location (cmd/ini)
     - qualifier
     - example
   * - database / d
     - Specify any number of FASTA formatted databases you would like to configure for EnTAP. Not necessary if you already have DIAMOND configured databases (.dmnd extension)
     - cmd
     - multi-string
     - /path/to/diamond/database
   * - threads / t
     - Specify thread number
     - cmd
     - number
     - 4
   * - out-dir
     - Specify the output directory for databases (EnTAP, DIAMOND, EggNOG) to be written to
     - cmd
     - string
     - /path/to/entap_output
   * - data-generate
     - Specify this flag if you would like to generate EnTAP databases rather than downloading from FTP (default). Only needed if you are having issues connecting to FTP.
     - cmd
     - flag
     - data-generate
   * - data-type
     - Specify which EnTAP database you'd like to generate/download. There are two provided in case there are any compatibility issues across platforms.
     
         * 0. Binary Database (default) - Recommended
         * 1. SQL Database - Slower although may resolve certain compatibilty issues across systems with low available memory
       This can be used multiple times to download/configure multiple databases at once. 
     - ini
     - multi-number
     - 0
 
