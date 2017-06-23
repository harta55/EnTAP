.. _rsem: https://github.com/deweylab/RSEM
.. _InterProScan: http://www.ebi.ac.uk/interpro/interproscan.html
.. _eggnog: https://github.com/jhcepas/eggnog-mapper
.. _diamond: https://github.com/bbuchfink/diamond
.. _GeneMarkS-T: http://exon.gatech.edu/GeneMark/

enTAP Introduction
==================

*enTAP* is a Eukaryotic Non-Model Transcriptome Annotation Pipeline designed to be fast, efficient, and user-friendly. It provides integration to several aspects of annotation and filtering analysis that will ultimately result in an more accurate match of sequences to gene ontology, pathway, and functional information.

Pipeline Stages:
----------------
#. Frame Selection
#. Expression Filtering
#. Similarity Search
#. Gene Ontology Integration

The software in the above stages are packaged within the enTAP repository (with the exception of GeneMarkS-T for frame selection). Installation and usage of enTAP can be seen in the sidebar.


Citations:
----------
[1] B. Buchfink, Xie C., D. Huson, "Fast and sensitive protein alignment using 
      DIAMOND", Nature Methods 12, 59-60 (2015).

[2] eggNOG 4.5: a hierarchical orthology framework with improved functional
      annotations for eukaryotic, prokaryotic and viral sequences. Jaime
      Huerta-Cepas, Damian Szklarczyk, Kristoffer Forslund, Helen Cook, Davide
      Heller, Mathias C. Walter, Thomas Rattei, Daniel R. Mende, Shinichi
      Sunagawa, Michael Kuhn, Lars Juhl Jensen, Christian von Mering, and Peer
      Bork. Nucl. Acids Res. (04 January 2016) 44 (D1): D286-D293. doi:
      10.1093/nar/gkv1248

[3] Fast genome-wide functional annotation through orthology assignment by
      eggNOG-mapper. Jaime Huerta-Cepas, Damian Szklarczyk, Lars Juhl Jensen,
      Christian von Mering and Peer Bork. Submitted (2016).

[4] Li, B., Dewey, C., & Liu, P. RSEM. In.

Software contained or used within this pipeline:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* `RSEM`_
* `DIAMOND`_
* `EggNog`_
* `InterProScan`_
* `GeneMarkS-T`_