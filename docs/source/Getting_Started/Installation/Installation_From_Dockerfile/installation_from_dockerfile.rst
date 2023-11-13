Installation from Dockerfile
==============================

If you'd prefer to use the Docker image to execute EnTAP, this can be done through :ref:`Docker<docker-label>` or :ref:`Singularity<singularity-label>`.

.. _docker-label:

Docker
----------------------

The Docker image can either be created locally from the Dockerfile, or pulled from Dockerhub. 

To download the latest image from Dockerhub:

.. code-block :: bash

    docker image pull plantgenomics/entap:latest
	
To download a specific version of the EnTAP Docker image:

.. code-block :: bash

	docker image pull plantgenomics/entap:vX.Y.Z

.. _singularity-label:

The list of arguments that can be ran through the image will be described in the following pages.

.. note:: Pay special attention to the entap_config.ini file located under /docker/entap_config.ini in the repo. The software execution paths here are specific to the Docker image paths and should not be changed!

Singularity
----------------------

To download the latest Docker image through Singularity:

.. code-block :: bash

    singularity pull entap.sif docker://plantgenomics/entap:latest

To download a specific version of the EnTAP Docker image through Singularity:

.. code-block :: bash

    singularity pull entap.sif docker://plantgenomics/entap:vX.Y.Z
	
The list of arguments that can be ran through the image will be described in the following pages.	

.. note:: Pay special attention to the entap_config.ini file located under /docker/entap_config.ini in the repo. The software execution paths here are specific to the Docker image paths and should not be changed!
