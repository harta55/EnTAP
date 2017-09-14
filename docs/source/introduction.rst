.. _rsem: https://github.com/deweylab/RSEM
.. _InterProScan: http://www.ebi.ac.uk/interpro/interproscan.html
.. _eggnog: https://github.com/jhcepas/eggnog-mapper
.. _diamond: https://github.com/bbuchfink/diamond
.. _GeneMarkS-T: http://exon.gatech.edu/GeneMark/

EnTAP Introduction
==================

The Eukaryotic Non-Model Transcriptome Annotation Pipeline (*EnTAP*) is designed to improve the accuracy, speed, and flexibility of functional gene annotation for de novo assembled transcriptomes in non-model eukaryotes. 

This software package addresses the fragmentation and related assembly issues that result in inflated transcript estimates and poor annotation rates.  Following filters applied through assessment of true expression and frame selection, open-source tools are leveraged to functionally annotate the translated proteins. 

Downstream features include fast similarity search across multiple databases, protein domain assignment, orthologous gene family assessment, Gene Ontology term assignment, and KEGG pathway annotation.  

The final annotation integrates across multiple databases and selects an optimal assignment from a combination of weighted metrics describing similarity search score, taxonomic relationship, and informativeness.  Researchers have the option to include additional filters to identify and remove potential contaminants and prepare the transcripts for enrichment analysis.  This fully featured pipeline is easy to install, configure, and runs much faster than comparable functional annotation packages.  It is developed to contend with many of the issues in existing software solutions.  

EnTAP is optimized to generate extensive functional information for the gene space of organisms with limited or poorly characterized genomic resources.


Pipeline Stages:
----------------
    * Transcriptome Filtering: designed to remove assembly artifacts and identify true CDS (complete and partial genes)
    1. Expression Filtering (RSEM)
    2. Frame Selection (GeneMARKS-T)

    * Annotation
    3. Similarity Search: optimized search against user-selected databases (DIAMOND).  
    4. Contaminant Filtering and Best Hit Selection: selects final annotation and identifies potential contaminants
    5. Orthologous Group Assignment: independent assignment of translated protein sequences to gene families (eggNOG).  Includes protein  domains (SMART/Pfam), Gene Ontology (GO) terms, and KEGG pathway assignment.

All of the software integrated into this pipeline are packaged within the EnTAP repository with the exception of GeneMarkS-T. Installation and usage of EnTAP is documented in this guide.

Citations:
----------
[1] A. Mitchell, H.-Y. Chang, and L. Daugherty, "The InterPro protein families database: the
      classification resource after 15 years," Nucleid Acids Research, vol. 43, no. D1, 2015.

[2] B. Buchfink, Xie C., D. Huson, "Fast and sensitive protein alignment using 
      DIAMOND", Nature Methods 12, 59-60 (2015).

[3] eggNOG 4.5: a hierarchical orthology framework with improved functional
      annotations for eukaryotic, prokaryotic and viral sequences. Jaime
      Huerta-Cepas, Damian Szklarczyk, Kristoffer Forslund, Helen Cook, Davide
      Heller, Mathias C. Walter, Thomas Rattei, Daniel R. Mende, Shinichi
      Sunagawa, Michael Kuhn, Lars Juhl Jensen, Christian von Mering, and Peer
      Bork. Nucl. Acids Res. (04 January 2016) 44 (D1): D286-D293. doi:
      10.1093/nar/gkv1248

[4] Fast genome-wide functional annotation through orthology assignment by
      eggNOG-mapper. Jaime Huerta-Cepas, Damian Szklarczyk, Lars Juhl Jensen,
      Christian von Mering and Peer Bork. Submitted (2016).

[5] Li, B., Dewey, C., & Liu, P. RSEM. In.

[6] P. Jones et al., "InterProScan 5: genome-scale protein function classification," (in eng),
      Bioinformatics, vol. 30, no. 9, pp. 1236-40, May 2014.

[7] S. Tang, A. Lomsadze, and M. Borodovsky, "Identification of protein coding regions in
      RNA transcripts," Nucleic Acids Research, 2015

Software contained or used within this pipeline:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* `RSEM`_
* `DIAMOND`_
* `EggNOG`_
* `GeneMarkS-T`_
* `InterProScan`_
