Frame Selection Flags
=============================

These are the flags specific to Frame Selection using either GeneMarkS-T or TransDecoder.

General Flags
------------------

*-*-frame-selection [integer] [INI]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Specify the frame selection software to use
    * 1. GeneMarkS-T 
    * 2. TransDecoder (default)

*-*-complete [INI]
^^^^^^^^^^^^^^^^^^^^^^
* Tell EnTAP to mark all of the transcripts as 'complete'. This will only be seen in the final output and will not affect the run.

TransDecoder Specific Flags
----------------------------------

*-*-transdecoder-m
^^^^^^^^^^^^^^^^^^^^^^^^
* Specify the minimum protein length
