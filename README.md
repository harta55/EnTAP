# EnTAP

![entap](docs/source/ENTAP_white_50.jpg?raw=true)

EnTAP (Eukaryotic Non-Model Transcriptome Annotation Pipeline) was designed to improve the accuracy, speed, and flexibility of functional gene annotation for de novo assembled transcriptomes in non-model eukaryotes.  This software package addresses the fragmentation and related assembly issues that result in inflated transcript estimates and poor annotation rates, while focusing primarily on protein-coding transcripts.  Following filters applied through assessment of true expression and frame selection, open-source tools are leveraged to functionally annotate the translated proteins.  

Downstream features include fast similarity search across three repositories, protein domain assignment, orthologous gene family assessment, and Gene Ontology term assignment.  The final annotation integrates across multiple databases and selects an optimal assignment from a combination of weighted metrics describing similarity search score, taxonomic relationship, and informativeness.  Researchers have the option to include additional filters to identify and remove contaminants, identify associated pathways, and prepare the transcripts for enrichment analysis.  

This fully featured pipeline is easy to install, configure, and runs significantly faster than comparable annotation packages.   It is developed to contend with many of the issues in existing software solutions.  EnTAP is optimized to generate extensive functional information for the gene space of organisms with limited or poorly characterized genomic resources.

Full Documentation can be found at:
http://entap.readthedocs.io/en/latest/

For information/bug reports, contact Alexander Hart at entap.dev@gmail.com

---
Copyright 2017-2023
Alexander Hart, Dr. Jill Wegrzyn, Dr. Stephen Ficklin, Josh Burns

EnTAP is protected under the GNU General Public License Version 3

Hart AJ, Ginzburg S, Xu M, et al. EnTAP: Bringing faster and smarter functional annotation to non-model eukaryotic transcriptomes. Mol Ecol Resour. 2020;20:591–604. https://doi.org/10.1111/1755-0998.13106
