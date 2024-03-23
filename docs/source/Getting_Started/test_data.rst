Test Data
==================
Before continuing to running EnTAP, it is advised to do a test run to ensure that everything is properly configured. There should be no errors in the test run. The test data resides within the |test_dir| directory of the main EnTAP directory. This will walk you through configuring a database for DIAMOND (if you haven't already done so) and executing EnTAP with and without frame selection. 

Before we begin, make sure that the paths in the |config_file_format| and |run_ini_file_format| files are correct. This test run will show you show to configure DIAMOND databases so we don't want to waste time and re-download the others (EggNOG and EnTAP) if they've already been downloaded. 

In |run_ini_file_format|, change the following lines:

.. code-block:: bash

    database=/test_data/uniprot_sprot.pep
    out-dir=test_data

Once the ini files have been setup, execute the following command to configure the test DIAMOND database:

.. code-block:: bash

    EnTAP --config --run-ini path/to/entap_config.ini --entap-ini path/to/entap_run.params


This should finish very shortly without any errors and you should find a uniprot_sprot.dmnd file within the |test_dir| directory. 

Next up is verifying the main execution stage! Once again, first ensure that the Config File has all of the correct paths. We are going to check an execution with and without frame selection. If you are not going to use frame selection, you may skip this test!

.. note:: The following tests will take longer as they will be testing the entire pipeline and running against the larger EggNOG database.

In |run_ini_file_format|, change the following lines for our frame selection test:

.. code-block:: bash
    input=/test_data/nucleotide_trinity.fnn
    database=/test_data/bin/uniprot_sprot.dmnd

To test EnTAP with the frame selection portion of the pipeline, execute the following command with the previous ini file changes:

.. code-block:: bash

    EnTAP --runP --run-ini path/to/entap_config.ini --entap-ini path/to/entap_run.params

Update the |run_ini_file_format| file again for our test without frame selection:

.. code-block:: bash

    input=/test_data/protein_transdecoder.faa
	
To test EnTAP without the frame selection portion of the pipeline, execute the following command with the previous ini file changes:

.. code-block:: bash

    EnTAP --runP --run-ini path/to/entap_config.ini --entap-ini path/to/entap_run.params

These should run without error and you should have several files within the created |out_dir| directory.

If any failures were seen during the above executions, be sure to go through each stage of installation and configuration to be sure everything was configured correctly before continuing. If you have no received any errors, the following pages will go through the main annotation portion of the pipeline.

.. note:: The |config_file| file may be distributed to each user to ensure they are pointing to the correct database and execution paths.
