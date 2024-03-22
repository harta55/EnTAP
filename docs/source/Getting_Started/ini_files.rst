EnTAP Ini Files
=====================
Before starting to Configure EnTAP, it's important to understand the ini files EnTAP uses and to setup the proper paths. EnTAP relies on several accompanying software packages and databases in order to run properly. This requires an easy way to specify parameters and paths rather than overloading the command line with a bunch of arguments. To resolve this, separate ini files are mandatory to be specified when running EnTAP as of version 1.1.0:

#. :ref:`|config_file|<entap-config-label>`
#. :ref:`|run_ini_file|<entap-params-label>`

.. _entap-config-label:

|config_file| File
---------------------------------------
The |config_file| contains database paths and execution paths for accompanying software. The intention is that this can be setup once and distributed to all users in a shared cluster environment, or all users of the EnTAP installation. It should likely not change between separate runs of EnTAP, unlike the |run_ini_file|. 

Many of the paths for accompanying software can simply use the globally installed command, such as :file:`diamond`. Or they full path to the executable. Below is the default |config_file| located in the EnTAP repository with it's defaults:

.. literalinclude:: ../../|config-file|

.. _entap-params-label:

|run_ini_file| File
---------------------------------------
The |run_ini_file| contains paths and parameters specific to each EnTAP run. The intention is that this parameter file will change between separate runs of EnTAP, either utilizing different databases or different parameters/cutoffs. Unlike the |config_file| this will be updated often by the user!

EnTAP gives a lot of power to the user to customize their run. However, defaults are provided for all, non-mandatory, parameters.  Below is an example of the |run_ini_file|:

.. literalinclude:: ../../|run_ini-file|

.. warning:: Be sure to at least set the DIAMOND and database flags before moving on as this is used during Configuration if you would like to configure databases for DIAMOND