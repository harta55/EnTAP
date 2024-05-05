.. |exp_dir| replace:: :file:`/expression/RSEM`
.. |exp_proc_dir| replace:: :file:`/expression/RSEM/processed`

Expression Analysis
=============================
The goal of expression filtering, or transcript quantification, is to determine the relative abundance levels of transcripts when taking into account the sequenced reads and how they map  back to the assembled transcriptome and using this information to filter out suspect expression profiles possibly originated from poor or incomplete assemblies. Filtering is done through the use of RSEM to determine the FPKM (fragments per kilobase per of million mapped reads) value. The FPKM , or a measurable number of expression, is then used as a threshold by EnTAP to filter the transcriptome. Anything lower than that threshold is removed.

Running Expression Analysis
----------------------------------
As mentioned above, RSEM (|rsem_github|) is utilized to determine the FPKM for each transcript. In order to run Expression Analysis, an alignment (BAM/SAM) file must be input to EnTAP through the :file:`align` flag in the |run_ini_file_format| file. A BAM file is preffered, but if a SAM file is input, RSEM/EnTAP will conver this over for use with EnTAP through :file:`convert-sam-for-rsem` RSEM package. Be sure to review the following relevant commands as well to tailor the FPKM threshold to what you prefer.

Expression Analysis Commands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. list-table:: **Expression Analysis Flags**
   :align: left
   :widths: 10 50 10 10 10 
   :header-rows: 1    
   
   * - param
     - description
     - location (cmd/R-ini,E-ini)
     - qualifier
     - example
   * - align / a
     - Path to the alignment file (either SAM or BAM format). Ignoring this flag will skip expression filtering. Be sure to look at the other Expression Analysis flags if using this.
     - R-ini
     - string
     - /path/to/alignment.bam
   * - rsem-calculate-expression
     - Specify the path to the :file:`rsem-calculate-expression` file generated during installation of RSEM
     - E-ini
     - string
     - /libs/RSEM-1.3.3/rsem-calculate-expression
   * - rsem-sam-validator
     - Specify the path to the :file:`rsem-sam-validator` file generated during installation of RSEM
     - E-ini
     - string
     - /libs/RSEM-1.3.3/rsem-sam-validator
   * - rsem-prepare-reference
     - Specify the path to the :file:`rsem-prepare-reference` file generated during installation of RSEM
     - E-ini
     - string
     - /libs/RSEM-1.3.3/rsem-prepare-reference
   * - rsem-convert-sam-for-rsem
     - Specify the path to the :file:`rsem-convert-sam-for-rsem` file generated during installation of RSEM
     - E-ini
     - string
     - /libs/RSEM-1.3.3/rsem-convert-sam-for-rsem
   * - fpkm
     - Specify the FPKM (fragments per kilobase of exon per million mapped fragments) cutoff for Expression Filtering. All transcripts below this number will be filtered out and removed.
     - R-ini
     - float
     - 0.5
   * - single-end
     - Signify your reads are single-end for RSEM execution instead of paired-end (default)
     - R-ini
     - bool
     - true
	 

Interpreting the Results
-----------------------------
The |exp_dir| folder will contain all of the relevant information for this stage of the pipeline. This includes many files generated from RSEM as well as files generated from EnTAP. Files generated from EnTAP are contained within the |exp_proc_dir| directory. 

RSEM generates many files, but the :file:`genes.results` file is what we are particularly interested in from the RSEM output. This contains the relevant FPKM values used for thresholding. The following files can be found within the |exp_dir| directory using an example input transcriptome titled "Species.fasta":

.. list-table:: **Expression Analysis Results**
   :align: left
   :widths: 10 50 10
   :header-rows: 1    
   
   * - filename
     - description
     - directory
   * - :file:`Species.fasta.genes.results`
     - Generated from RSEM. Contains important information for each transcsript such as FPKM and TPM
     - |exp_dir|
   * - :file:`Species_removed.fasta`
     - Generated from EnTAP. Contains all of the transcripts that have been filtered out due to having an FPKM value below the threshold
     - |exp_proc_dir|
   * - :file:`Species_kept.fasta`
     - Generated from EnTAP. Contains all of the transcripts that have been retained due to having an FPKM threshold above the user input one
     - |exp_proc_dir|


Expression Analysis Headers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TSV files generated from EnTAP will have the following headers from Expression Analysis.

    * FPKM
    * TPM
    * Expression Effective Length
