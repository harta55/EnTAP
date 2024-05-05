.. _configuration_flags-label:

Configuration Flags
============================

These are the flags for the configuration process of EnTAP. These will be used via the command line (denoted CMD), |run_ini_file_format| file (denoted R-ini), or |config_file_format| file (denoted E-ini).

There are a few data types (qualifiers) to keep in mind used throughout these ini files. Anything that is specifies as a 'multi' type means that the parameter may be entered multiple times. If it is in an ini file, each parameters must be separated by a comma (','). Example for multi-integer:"1,2,3" (entered without quotes)

.. list-table:: **Required/Recommended Flags**
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1
   
   * - param
     - description
     - location (cmd/R-ini,E-ini)
     - qualifier
     - example
   * - config
     - Must be specified in order for EnTAP to run Configuration as opposed to normal execution
     - cmd
     - flag
     - config
   * - run-ini
     - Point to the |run_ini_file_format| to specify database paths to be configured for DIAMOND, thread usage, and more.
     - cmd
     - string
     - /path/to/|run_ini_file|
   * - entap-ini
     - Point to the |config_file_format| to specify database paths to prevent redownloading and DIAMOND executable in order to format FASTA databases for DIAMOND.
     - cmd
     - string
     - /path/to/|config_file|
   * - eggnog-map-data
     - Specify the output directory for the EggNOG database (|eggnog_map_sql_db_file_format|). EnTAP will check for the eggnog.db database within this specified directory. If it is not there, it will be redownloaded. So set this to prevent duplicate downloads!
     - E-ini
     - string
     - /path/to/eggnog_database_directory
   * - eggnog-map-dmnd
     - Specify the filepath to the eggnog_proteins.dmnd EggNOG database. EnTAP will check this path for the database. If it does not exist, it will be downloaded. So set this to prevent duplicate downloads!
     - E-ini
     - string
     - /path/to/eggnog_proteins.dmnd
   * - entap-db-bin
     - Path to the |entap_db_bin_file_format| database. Either this or the EnTAP SQL database must be used. The binary database is the default and is recommended. EnTAP will check for this path during Configuration. If it is not there, it will be redownloaded. So set this to prevent duplicate downloads!
     - E-ini
     - string
     - /path/to/|entap_db_bin_file|
   * - diamond-exe
     - Specify the execution method for DIAMOND. This can be a path to the :file:`diamond` file generated during installation, or simply the command if installed globally. Used to configure databases for EnTAP to use during Similarity Searching
     - E-ini
     - string
     - diamond
   * - database / d (optional)
     - Specify any number of FASTA formatted databases you would like to configure for EnTAP. Not necessary if you already have DIAMOND configured databases (.dmnd extension)
     - R-ini
     - multi-string
     - /path/to/diamond/database
   * - threads / t
     - Specify thread number
     - R-ini
     - number
     - 4
   * - out-dir
     - Specify the output directory for databases (EnTAP, DIAMOND, EggNOG) to be written to
     - R-ini
     - string
     - /path/to/entap_output

     
.. list-table:: **Optional Flags**
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1
   
   * - param
     - description
     - location (cmd/R-ini,E-ini)
     - qualifier
     - example
   * - data-generate
     - Specify this flag if you would like to generate EnTAP databases rather than downloading from FTP (default). Only needed if you are having issues connecting to FTP.
     - R-ini
     - bool
     - true
   * - data-type
     - Specify which EnTAP database you'd like to generate/download. There are two provided in case there are any compatibility issues across platforms.
     
         * 0. Binary Database (default) - Recommended
         * 1. SQL Database - Slower although may resolve certain compatibilty issues across systems with low available memory
		 
       This can be used multiple times to download/configure multiple databases at once. 
     - R-ini
     - multi-number
     - 0
   * - entap-db-sql
     - Path to the |entap_db_sql_file_format| database. Either this or the EnTAP binary database must be used. EnTAP will check for this path during Configuration. If it is not there, it will be redownloaded. So set this to prevent duplicate downloads!
     - E-ini
     - string
     - /path/to/|entap_db_sql_file|