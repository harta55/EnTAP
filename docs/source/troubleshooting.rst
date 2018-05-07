Troubleshooting
==================
This page contains several useful ways to troubleshoot whenever you run into an issue with EnTAP execution or installation. If you do not come to a solution through checking inputs, paths, or the configuration file I'd review these in the order I am laying them out.

#. :ref:`Error Codes<error_codes-label>` : Provide specific information on where your run failed within the EnTAP pipeline. This could be during execution or installation. These can be as general as 'DIAMOND has failed,' or much more detailed. They will be more detailed if the failure was detected by EnTAP rather than accompanying software.

#. :ref:`.err and .out Files<err_out-label>` : Provide information from the software that has failed. These are essentially what the software would print out to the user if an error occurs. Thus, they are out of EnTAP's control and are the most useful item when troubleshooting why a specific program failed. This should be the second thing you look for after checking the standard error code produced by EnTAP. An example would look like: 'DIAMOND was unable to finish execution due to a memory overflow issue.' 

#. :ref:`FAQs<faq-label>` : If you're still running into trouble, check out these common issues!


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

.. _faq-label:

FAQs
-------------------

#. I'm having issues reading the EnTAP database!

        * First, ensure you are specifying the correct paths in the configuration file. Next, attempt to generate the database locally (through - - data-generate flag during configuration). This is likely the result of incompatibility in Boost versions. The SQL version should always be compatible with your system as well, although slower.

#. Eggnog-mapper is failing when running DIAMOND.

        * Eggnog-mapper leverages DIAMOND to search against the EggNOG databases. In order to do this, it uses a global call to DIAMOND. This may not work if you do not have it installed globally. Additionally, the DIAMOND EggNOG database may not be compatible with your local DIAMOND version. If so, you'll have to re-index the fasta version of this database (found from the EggNOG FTP) with your version of DIAMOND. 