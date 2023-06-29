API Flags
=====================

These flags are used as part of an ongoing effort to create an EnTAP API to support software packages incorporating EnTAP within their framework. 

*-*-api-taxon [string] [CMD]
------------------------------
* Check whether a species is part of the EnTAP Taxonomy database specified in the ini file
* Returns json formatted text indicating whether the species was found
* Format **must** replace all spaces with underscores ('_') as follows: "- -taxon homo_sapiens" or "- -taxon primates"