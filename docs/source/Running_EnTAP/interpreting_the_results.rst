.. |final_dir| replace:: */final_results*
.. |transc_dir| replace:: */transcriptomes*


.. _GO: http://www.geneontology.org/


Interpreting the Results
=================================

*EnTAP* provides many output files at each stage of execution to better see how the data is being managed throughout the pipeline:

#. :ref:`Final Annotation Results<final-label>`
#. :ref:`Log File / Statistics<log-label>`
#. :ref:`Transcriptomes<transc-label>`
#. :ref:`Expression Filtering
#. :ref:`Frame Selection
#. :ref:`Similarity Searching
#. :ref:`Orthologous Groups/Ontology
#. :ref:`Protein Families (optional)

The two files to check out first are the :ref:`final annotations<final-label>` and :ref:`log file<log-label>`. These files contain a summary of all the information collected at each stage, including statistical analyses. The remaining files are there for a more in depth look at each stage. All files will be contained in "entap_outfiles" directory as default, or different if the - - out-dir flag was specified.

.. _final-label:

Final Annotations
-----------------------

The final EnTAP annotations are contained within the |final_dir| directory. These files are the summation of each stage of the pipeline and contain the combined information. So these can be considered the most important files! The "full_entap.tsv" file will contain all of the information gathered throughout the pipeline summarized in one file. This will include annotated, unannotated, and contaminated sequences. 

All .tsv files in this section may have the following header information (from left to right) separated by each portion of the pipeline. Some headers will not be shown if that part of the pipeline was skipped or the information was not found for any of the input sequences. TSV formatted files support Tidyverse format (including 'NA' being used for empty data cells).

General Header Information
    * Query sequence ID

Frame Selection Header Information (optional)
    * Open Reading Frame

Expression Analysis Header Information (optional)
    * FPKM
    * TPM
    * Effective Length

Similarity Search Header Information
    * Subject sequence ID
    * Percentage of identical matches
    * Alignment length
    * Number of mismatches
    * Number of gap openings
    * Start of alignment in query
    * End of alignment in query
    * Start of alignment in subject
    * End of alignment in subject
    * Expect (e) value
    * Query coverage
    * Subject title
    * Species
    * Taxonomic Lineage
    * Origin Database
    * Contaminant (yes/no if the hit was flagged as a contaminant)
    * Informative (yes/no if he hit was flagged as informative)

Similarity Search UniProt Header Information (optional if aligning against SwissProt database)
    * UniProt Database Cross References
    * UniProt Additional Information
    * UniProt KEGG Terms
    * UniProt GO Biological
    * UniProt GO Cellular
    * UniProt GO Molecular

Ontology EggNOG Header Information
    * Seed Ortholog
    * Seed E-Value
    * Seed Score
    * Predicted Gene
    * Taxonomic Scope
    * OGs (orthologous groups assigned)
    * EggNOG Description (EggNOG)
    * KEGG Terms (EggNOG)
    * GO Biological (Gene Ontology)
    * GO Cellular (Gene Ontology)
    * GO Molecular (Gene Ontology)
    * BIGG Reaction

Ontology InterProScan Header Information
    * IPScan GO Biological 
    * IPScan GO Cellular
    * IPScan GO Molecular
    * Pathways
    * InterPro (InterPro database entry)
    * Protein Database (database assigned. Ex: pfam)
    * Protein Description (description of database entry)
    * E Value (E-value of hit against protein database)


    * full_entap.tsv

        * This .tsv file is essentially a final report from EnTAP that will have the headers as mentioned previously, summarizing the results of the entire pipeline
        * Since this includes every single transcript, there will be annotated, unannotated, and contaminated sequences. Further filtering of transcripts (for example if you are only interested in those transcripts that were annotated) can be done with this file or the below files

    * annotated.faa / .fnn / .tsv

        * Nucleotide/protein fasta files along with tsv file containing all sequences that either align databases through similarity searching or through the ontology stage

    * unannotated.faa / .fnn / .tsv

        * Nucleotide/protein fasta files along with tsv file containing all sequences that did not align either through similarity searching nor through the ontology stage

    * annotated_contam.faa / .fnn / .tsv

        * Nucleotide/protein fasta files along with tsv file containing all annotated sequences that were flagged as a contaminant

    * annotated_without_contam.faa / .fnn / .tsv

        * Nucleotide/protein fasta files along with tsv file containing all annotated sequences that were not flagged as a contaminant

    * x_enrich_geneid_go.tsv

        * Tab-deliminated file that can be used for Gene Enrichment
        * First column contains the gene ID and second column contains the Gene Ontology term corresponding to the gene ID

    * x_enrich_geneid_len.tsv

        * Tab-deliminated file that can be used for Gene Enrichment
        * First column contains the gene ID and second columns contains the effective length from Expression Analysis. This file will not be printed if Expression Analysis has not been ran
        * Note: the Length column will not be printed when Expression Filtering has not been performed

    * x_gene_ontology_terms.tsv

        * Tab-deliminated file that can be used for Gene Enrichment
        * Columns are as follows: Sequence ID, Gene Ontology Term ID, Gene Ontology Term, Gene Ontology Category, and Effective Length
        * Note: the Effective Length column will not be printed when Expression Filtering has not been performed
		

.. _log-label:

Log File / Statistics
-----------------------------

The log file contains a statistical analysis of each stage of the pipeline that you ran. I'll give a brief outline of some of the stats performed:

#. Initial Statistics

    * Transcriptome statistics: n50, n90, average gene length, longest/shortest gene
    * Summary of user flags
    * Summary of execution paths (from config file)

#. Expression analysis

    * Transcriptome statistics: n50, n90, average gene length, longest/shortest gene
    * Summary of sequences kept/removed after filtering

#. Frame Selection

    * Transcriptome statistics: n50, n90, average gene length, longest/shortest gene
    * Summary of frame selection: Partial, internal, complete genes. Genes where no frame was found

#. Similarity Searching

    * Contaminant/uninformative/informative count
    * Phylogenetic/contaminant distribution of alignments
    * Alignment distribution based upon frame results (partial/internal/complete)
    * Sequence count that did not align against a database reference
    * Statistics calculated for each individual database and final results

#. Gene Family Assignment

    * Phylogenetic distribution of gene family assignments
    * Gene Ontology category distribution (biological processes, molecular function, cellular component)

#. InterProScan

    * Additional statistics coming soon!

#. Final Annotation Statistics

    * Statistical summary of each stage
    * Runtime


.. _transc-label:

Transcriptomes
---------------------
The |transc_dir| contains the original, processed, and final transcriptomes being used by EnTAP. The files are as follows with the 'transcriptome' tag based upon the name of your input transcriptome:

* transcriptome.fasta

    * This file is essentially a copy of your input transcriptome. The sequence ID's may be changed depending on whether you selected the 'trim' flag or otherwise.

* transcriptome_expression_filtered.fasta

    * As the name implies, this transcriptome is the resultant of the Expression Filtering stage with sequences removed that fall under the FPKM threshold you have specified.

* transcriptome_frame_selected.fasta

    * This transcriptome is the resultant of Frame Selection. Sequences in which a frame was not selected are removed and those with a frame are kept in this file. As a result, this file will always be in protein format. 

* transcriptome_final.fasta

    * This is your final transcriptome following the "Transcriptome Filtering" stage of EnTAP. **This transcriptome will be used for the later stages of the pipeline** (Similarity Searching and Ontology). Depending on which methods of execution you chose (runN / runP), the result here may be either protein or nucleotide with Frame Selection and/or Expression Filtering.