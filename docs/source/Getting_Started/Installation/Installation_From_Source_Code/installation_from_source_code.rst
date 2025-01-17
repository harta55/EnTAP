.. _Boost: http://www.boost.org/users/download/
.. _Perl: https://www.perl.org/
.. _Python: https://www.python.org/
.. _RSEM: https://github.com/deweylab/RSEM
.. _EggNOG-mapper: https://github.com/jhcepas/eggnog-mapper
.. _DIAMOND: https://github.com/bbuchfink/diamond
.. _CMake: https://cmake.org/
.. _InterProScan: https://github.com/ebi-pf-team/interproscan
.. _TransDecoder: https://github.com/TransDecoder/TransDecoder/releases
.. _NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/taxonomy

Installation from Source Code
=================================

#. :ref:`System Requirements<sys-label>`
#. :ref:`Download EnTAP Source Code<ent-label>`
#. :ref:`Installing Pipeline Software<dep-label>`
#. :ref:`EnTAP Installation<entap-label>`

.. _sys-label:

System Requirements
-----------------------------------
  
    * Operating System

        * UNIX-based systems
        * Tested on 64 bit systems: ubuntu 16.04, Rocks 6.1, Centos 6.3

    * Storage Minimum

        * EnTAP Database (Gene Ontology References + UniProt Mapping + NCBI Taxonomy): 1.5Gb
        * EggNOG Databases: 24Gb
        * DIAMOND Databases: ~13Gb (with RefSeq Complete Protein + Uniprot Swiss-Prot)
        * Additional storage for files generated depending on transcriptome size: upwards of 15Gb

    * Memory

        * At least 16 Gb of RAM (will vary depending on DIAMOND database sizes). More memory is highly recommended to reduce execution times.

.. _dep-label:

Dependency Check
-----------------------------------
Before continuing on in the installation process, ensure that the following dependencies are fully installed on your system:

    * C++11 compiler (GCC 4.8.1 or later)
	
    * CMake_ (3.00 or later)
	
    * Python_ (2.7.12 or later) with support for the following modules	
    
        * Matplotlib (figures generated by EnTAP)
		
    * Unix wget (generally included in most distros)
	
    * Unix gzip/tar (generally included in most distros)


.. _ent-label:

Downloading EnTAP Source Code
----------------------------------------
First, download and extract the latest release(tagged) version from GitLab:
https://gitlab.com/EnTAP/EnTAP/tags

This repository contains the source code for RSEM, TransDecoder, and DIAMOND to be used with EnTAP.

.. _pipe-label:

Installing Pipeline Software
--------------------------------------------
EnTAP leverages several software distributions within the pipeline to provide the best quality annotations. The packages used (and their current/compatible versions) can be seen below. This is not to say that newer versions will not be compatible, however they have not been tested yet with EnTAP.

.. note:: If the software is already installed on your system, this stage can be skipped

Supported Software:
    * RSEM_ (Expression Filtering with alignment file): version 1.3.3 packaged with EnTAP

        * Version 1.3.0
        * Version 1.3.3

    * TransDecoder_ (Frame Selection): version 5.7.1 packaged with EnTAP
	
        * Version 5.3.0
        * Version 5.7.1

    * DIAMOND_ (Similarity Search): version 2.1.8 packaged with EnTAP

        * Version 2.0.10 (minimum required DIAMOND version)
        * Version 2.1.8
		
    * EggNOG-mapper_ (Gene Family assignment): version 2.1.12 packaged with EnTAP
	   
	    * Version 2.1.12

    * InterProScan_ (Protein Databases): must be installed separately
   
        * Version 5.19

.. _diamond-label:

DIAMOND Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^
From root EnTAP directory...

.. code-block :: bash

    cd libs/
    tar -xvzf diamond-v2.1.8.tar.gz
    cd diamond-v2.1.8
    mkdir bin
    cd bin
    cmake ..

Run the following command to install globally:

.. code-block :: bash

    make install

Run the following command to compile:

.. code-block :: bash

    make
	
All set! Ensure that DIAMOND has been properly setup and add the correct path to the |config_file| file. If installed globally, add 'diamond' (without quotes) to the file. If installed locally, add 'path/to/EnTAP/libs/diamond-2.1.8/bin/diamond'.
	
.. _eggnog-mapper-label:

EggNOG-mapper Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
EggNOG-mapper comes packaged with EnTAP, but can also be downloaded from the GitHub. If installing through the packaged version...

.. code-block :: bash

    cd libs/
    tar -xvzf eggnog-mapper-2.1.12.tar.gz
    cd eggnog-mapper-2.1.12

Run the following command to install globally:

.. code-block :: bash

    python setup.py install


All set! Ensure that EggNOG-mapper has been properly setup and add the correct path to the |config_file| file. If installed globally, add 'emapper.py' (without quotes) to the file. If installed locally, add 'path/to/EnTAP/libs/eggnog-mapper-2.1.12/emapper.py'.

.. _rsem-label:

RSEM Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^

From root EnTAP directory...

.. code-block :: bash

    cd libs/
    tar -xvzf RSEM-v1.3.3.tar.gz
    cd RSEM-v1.3.3
    make
    make ebseq

Run the following command to install globally:

.. code-block :: bash

    make install

All set! Ensure that RSEM has been properly setup and add the correct path to the entap_config.txt file. If installed globally keep blank. If installed locally, add 'path/to/EnTAP/libs/RSEM-1.3.0/'.

.. _entap-label:

EnTAP Installation
----------------------------

Once dependencies and pipeline software have been installed, you can now continue to install EnTAP! 

Within the main directory, execute the following command:

.. code-block :: bash

    cmake CMakeLists.txt

This will generate a MakeFile. Then execute:

.. code-block :: bash

    make

Or to install to a destination directory:

.. code-block :: bash

    cmake CMakeLists.txt -DCMAKE_INSTALL_PREFIX=/destination/dir

.. code-block :: bash

    make install

If you receive no errors, please move on to the last stage in installation, configuration.
