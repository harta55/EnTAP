Frame Selection Flags
=============================

These are the flags specific to Frame Selection using either GeneMarkS-T or TransDecoder. These will be used via the command line (denoted CMD) or ini file (denoted INI).

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

*-*-transdecoder-m [INI]
^^^^^^^^^^^^^^^^^^^^^^^^
* Specify the minimum protein length

*-*-transdecoder-no-refine-starts [INI]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Specify this flag if you would like to pipe the command, '--no_refine_starts' to TransDecoder. 
* "This will 'start refinement identifies potential start codons for 5' partial ORFs using a PWM, process on by default". (From TransDecoder documentation)
* Default: False
