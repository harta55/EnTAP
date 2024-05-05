.. |diamond_git| replace:: https://github.com/bbuchfink/diamond
.. |hgt_dir| replace:: :file:`/horizontal_gene_transfer/HGT_DIAMOND`
.. |hgt_proc_database_dir| replace:: :file:`/horizontal_gene_transfer/HGT_DIAMOND/processed/ref_database`
.. |hgt_proc_dir| replace:: :file:`/horizontal_gene_transfer/HGT_DIAMOND/processed/`

Horizontal Gene Transfer
===============================
Horizontal gene transfer (HGT) is an important evolutionary mechanism whereby organisms acquire foreign genetic material from other distantly related species. HGT detection can reveal key insights into how organisms adapt to new environments and acquire new functional capabilities.

EnTAP now provides the capability to predict horizontally acquired genes by comparing a genome against specified databases of potential "donor" and "recipient" species. This unlocks new opportunities for studying the impacts of HGT and how it enables evolutionary innovations. For purposes of HGT analysis, donor is the catch-all term for all of the databases taken into consideration that are likely to have been the origin of the horizontally transferred gene. Recipient databases encompass all potential (“closely” related) organisms that may have been the recipient of the horizontally transferred gene.

The HGT detection methodology utilizes sequence similarity searches (DIAMOND, |diamond_git|) to identify genes in the input genome with closer relationships to possible donor groups rather than the native genomic background. The detections can be further verified using genomic context and manual review.

HGTs identified in this tool are **unique HGTs** in the genome. Genes transferred between closely related species are also considered HGT (i.e., a gene transferred between Streptophyta and Spermatophyta), but the current version of EnTAP does not identify these cases.

.. note::
    An important pre-requisite to HGT analysis is having access to a reference genome with the annotation GFF file

Running Horizontal Gene Transfer Analysis
------------------------------------------------
The horizontal gene transfer (HGT) detection workflow consists of the following key steps:

    #. Configure DIAMOND databases - The user specifies "donor" and "recipient" DIAMOND databases that represent potential sources and recipients of HGT respectively

    #. Run DIAMOND - The input genome is aligned against the donor and recipient databases using DIAMOND to identify hits
	
    #. :ref:`Identify HGT candidates<verify_hgt-label>` - Based on the DIAMOND alignments, genes are flagged as potential HGT events if they have hits to donor species but not recipient species
 
    #. :ref:`Verify HGT candidates<verify_hgt-label>` - The genomic context (neighboring genes) of each HGT candidate is extracted from the input GFF file. Their taxonomic classifications are checked to determine if the HGT is likely real or an artifact
	
To perform HGT analysis, the user must input several required commands. An example from the |run_ini_file_format| file can be seen below:

.. code-block:: bash
    
    hgt-recipient=/path/to/recipient_database1.dmnd,/path/to/recipient_database2.dmnd
    hgt-donor=/path/to/donor_database1.dmnd,/path/to/donor_database2.dmnd
    hgt-gff=/path/to/gff_example.gff
	
As an example for determining your donor/recipient databases - *Physcmitrellopsis africana* is a moss endemic to the coastal forests of South Africa. We considered the recipient databases to encompass all organisms under Streptophyta, Spermatophyta, and Tracheophyta. Donor databases on the other hand include all organisms under archaea, bacteria, metazoan, and viral databases. The selection of databases is dependent on your genome.
	
The GFF file must follow several requirements:
    * Protein identifiers must match between FASTA and GFF attribute fields
    * Primary transcripts only (longest isoform for each gene)
    * Feature type = 'transcript' or 'mRNA'
    * Must be in relative order. This can be accomplished if it is ran through software such as agat_sp_keep_longest_isoform (https://agat.readthedocs.io/en/latest/tools/agat_sp_keep_longest_isoform.html)

Horizontal Gene Transfer Commands
----------------------------------------------
.. list-table:: **Horizontal Gene Transfer Flags**
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1    
   
   * - param
     - description
     - location (cmd/R-ini,E-ini)
     - qualifier
     - example
   * - hgt-donor
     - Specify the DIAMOND configured (.dmnd extension) donor databases for Horizontal Gene Transfer analysis. Separate databases with a comma (',')
     - R-ini
     - multi-string
     - path/to/donor/database1.dmnd,path/to/donor/database2.dmnd
   * - hgt-recipient
     - Specify the DIAMOND configured (.dmnd extension) recipient databases for Horizontal Gene Transfer analysis. Separate databases with a comma (',')
     - R-ini
     - multi-string
     - path/to/recipient/database1.dmnd,path/to/recipient/database2.dmnd
   * - hgt-gff
     - Specify path to the GFF file for HGT analysis. The input GFF must satisfy the following:
	 
           * Protein identifiers must match between FASTA and GFF attribute fields
           * Primary transcripts only (longest isoform for each gene)
           * Feature type = 'transcript' or 'mRNA'
           * Must be in relative order. This can be accomplished if it is ran through software such as agat_sp_keep_longest_isoform (https://agat.readthedocs.io/en/latest/tools/agat_sp_keep_longest_isoform.html)
		   
     - R-ini
     - string
     - path/to/gff/file.gff
   * - diamond-exe
     - Specify the execution method for DIAMOND. This can be a path to the :file:`diamond` file generated during installation, or simply the command if installed globally. DIAMOND is leverage for HGT analysis
     - E-ini
     - string
     - diamond

.. _verify_hgt-label:

Identifying and Verifying HGT Candidates
----------------------------------------------
After DIAMOND has been ran against the donor and recipient databases, EnTAP will analyze this information to identify then verify HGT candidates. 

Identifying HGT Candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
HGT candidates are identified through the following process:

    #. If alignments against Donor databases > 0, but not all

    #. If no alignments against Recipient databases

If all steps above pass, a potential HGT candidate has been found. These genes are then verified.

Verifying HGT Candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
After obtaining a list of putative HGT genes based on homology searches, further verification is necessary to check if they are true events or false positives. This involves inspecting the genomic context of each candidate using the input GFF file. The following logis is employed:

    #. If neighbors exist for our HGT Candidate

    #. If neighboring genes are not HGT Candidates
	   
    #. If neighboring genes have no Donor database alignments?

If all of the above are true, we have a verified horizontally transferred gene!

Interpreting the Results
-----------------------------
The |hgt_dir| folder will contain all of the relevant information for this stage of the pipeline. This includes many files generated by DIAMOND as well as EnTAP. This directory highly resembles that of the directory within Similarity Searching with similar EnTAP files generated (|hgt_proc_dir|).

Below are example files with an input transcriptome labelled 'species' and a reference database labelled 'ref_database'.

.. list-table:: **Horizontal Gene Transfer Results**
   :align: left
   :widths: 10 50 10
   :header-rows: 1    
   
   * - filename
     - description
     - directory
   * - :file:`blastp_species_ref_database.out`
     - Generated from DIAMOND. Contains a lot of information from the DIAMOND search including e-value, coverage, reference database descriptions, and much more. This is a typical output file from a BLAST type of search. The filename is prepended with either blastp or blastx depending on if runP or runN was used. A file like this will be generated for each donor/recipient database used.
     - |hgt_dir|
   * - :file:`blastp_species_ref_database_std.err/.out`
     - Generated from DIAMOND. These files are will contain any error or general information produced from the DIAMOND run. A file like this will be generated for each donor/recipient database used.
     - |hgt_dir|
   * - :file:`hgt_diamond_annotated.faa/.fnn/.tsv`
     - Generated from EnTAP. Contains all of the best alignments (protein/nucleotide format) that were selected from this database based on e-value/coverage. Note: Protein or nucleotide information may not be available to report depending on your type of input sequences or runN vs. runP.
     - |hgt_proc_database_dir|
   * - :file:`unannotated.faa/.fnn/.tsv`
     - Generated from EnTAP. Sequences (protein/nucleotide) from the transcriptome that did not hit against this particular reference database. This does not include sequences that were lost during expression filtering or frame selection. Again this will be seen per database.
     - |hgt_proc_database_dir|
   * - :file:`hgt_diamond_unselected_hits.faa/.fnn/.tsv`
     - Generated from EnTAP. DIAMOND alignment can result in several hits for each query sequence. With only one best alignment being selected, the rest are unselected and end up here. Unselected hits can be due to a low e-value or coverage.
     - |hgt_proc_database_dir|
   * - :file:`hgt_candidates.faa/.fnn`
     - Generated from EnTAP. Contains FASTA formatted confirmed unique HGTs. 
     - |hgt_proc_dir|
   * - :file:`hgt_candidates.tsv`
     - Generated from EnTAP. Contains alls confirmed unique HGTs in TSV format with headers: Sequence ID | Donor Database | Species | Database Description. A new line will be generated for each donor database that each HGT aligned against. So there may be multiple lines with the same HGT! The *Species* and *Database Description* is taken from the donor database provided.   
     - |hgt_proc_dir|


Horizontal Gene Transfer Headers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TSV files generated from EnTAP will have the following headers from Horizontal Gene Transfer analysis.

    * Horizontally Transferred Gene
