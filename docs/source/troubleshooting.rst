Troubleshooting
==================
This page contains several useful ways to troubleshoot whenever you run into an issue with EnTAP execution or installation. If you do not come to a solution through checking inputs, paths, or the configuration file I'd review these in the order I am laying them out.

#. :ref:`Debug File<debug_file-label>` : EnTAP will pipe any errors to the Debug File printed after every execution as well as the 'standard' output. Depending on your system, the 'standard' output may be placed into separate files or printed to the terminal

#. :ref:`Error Codes<error_codes-label>` : Provide specific information on where your run failed within the EnTAP pipeline. This could be during execution or installation. These can be as general as 'DIAMOND has failed,' or much more detailed. They will be more detailed if the failure was detected by EnTAP rather than accompanying software.

#. :ref:`.err and .out Files<err_out-label>` : Provide information from the software that has failed. These are essentially what the software utilized by EnTAP (ex: DIAMOND) would print out to the user if an error occurs. Thus, they are out of EnTAP's control and may provide some further insight into the issue. 

#. :ref:`FAQs<faq-label>` : If you're still running into trouble, check out these common issues!

.. _debug_file-label:

Debug File
------------------------

The Debug File contains detailed information of the EnTAP run. The more relevant information will reside at towards the bottom of the file. This will explicitly tell you what went wrong with the execution. As mentioned above, the pertinent error information will be printed to the terminal as well as to this file. 

.. _error_codes-label:

Error Codes
------------------------

EnTAP has quite a few error codes that are generally going to have messages attached to them when your execution fails. I'll list the error codes and possible messages/troubleshooting ideas here. This is a continually evolving section where additional codes will be added as they are created (and I find the time to update). However, pre-existing codes will NOT be changed!

0. No error!

10. There is some error in your inputs to EnTAP. This is most likely caused by your command arguments. Make sure to check these out as well as any relevant paths. EnTAP should provide the specific error attached to this code as to why the run failed.

12. Some error in parsing the configuration file attached to EnTAP. Make sure it is following the proper format as when it was first downloaded from the GitLab page. There should be no spaces or extra line breaks after the '=' for each parameter.

13. This error is called when the configuration file could not be generated. If EnTAP does not recognize a configuration file, it will attempt to generate one for you. 

14. This error is called when a configuration file could not be found in either your input, or the execution directory for EnTAP. It will merely generate the file and exit execution to allow you to fill in any parameters you may want to include. 

15. There was an issue in download the EnTAP Taxonomic Database. This will only be called when you are running the - -config command. Some solutions can be: check your internet connection, ensure you have the required Python libraries installed (it currently uses Python as a means of download), and try running this outside of EnTAP (with the commmand python tax_download.py -o output_file.txt). If you decide to download the text version separately, just be sure to include that in the entap_config.txt file and it will convert it to a binary when you run - - config again.

.. _err_out-label:

.err and .out Files
---------------------

The .err and .out files contain all of the standard output from software utilized by EnTAP. Consider DIAMOND...EnTAP will utilize DIAMOND for Similarity Searching. This software has it's own set of error and output information. EnTAP will print these to the relevant .err and .out files within the proper directory (in this case similarity_search/DIAMOND). Again, this information will also be printed to the Debug File. 

.. _faq-label:

FAQs
-------------------

This list will be continually updated with common problems experienced by users.


------

The following issues remain for historical purposes. They are not relevant beyond EnTAP v0.9.0.

#. I'm having issues reading the EnTAP database!

        * First, ensure you are specifying the correct paths in the configuration file. Next, attempt to generate the database locally (through - - data-generate flag during configuration). This is likely the result of incompatibility in Boost versions. The SQL version should always be compatible with your system as well, although slower.

#. Eggnog-mapper is failing when running DIAMOND.

        * Eggnog-mapper leverages DIAMOND to search against the EggNOG databases. In order to do this, it uses a global call to DIAMOND. This may not work if you do not have it installed globally. Additionally, the DIAMOND EggNOG database may not be compatible with your local DIAMOND version. If so, you'll have to re-index the fasta version of this database (found from the EggNOG FTP) with your version of DIAMOND. 

#. Eggnog-mapper is failing during parsing of the data due to "too_many_columns"
         * This is due to an incompatibility between EnTAP and the version of Eggnog-mapper you are using. Please use the version provided in the EnTAP repository (0.7.4.1-beta). Additionally, replace the eggnog.db file you are using with the database at the following address: http://eggnogdb.embl.de/download/emapperdb-4.5.0/eggnog.db.gz. This incompability will be resolved shortly (if not already)