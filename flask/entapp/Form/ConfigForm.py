"""
Contains the ConfigForm class.
"""
from flask_wtf import FlaskForm
from wtforms import BooleanField
from wtforms import SelectField
from wtforms import SelectMultipleField
from wtforms import StringField
from wtforms import SubmitField
from wtforms.fields.html5 import DecimalField
from wtforms.fields.html5 import IntegerField




class ConfigForm(FlaskForm):
    """
    Detailed description.
    """
    dataType = SelectField(
        "Data Type"
        ,choices=[("0","Serialized Database"),("1","SQLite Database")]
        ,description="Specifies which database EnTAP will download/generate."
    )
    fpkm = DecimalField(
        "FPKM Threshold"
        ,description="Specifies the FPKM threshold with expression analysis. EnTAP will filter out "
                     "transcripts below this value!"
    )
    singleEnd = BooleanField("BAM/SAM files generated through single-end reads.")
    seqComplete = BooleanField("All sequences are complete proteins.")
    minProteinLength = IntegerField(
        "Minimum Protein Length"
        ,description="Specifies the minimum protein length."
    )
    refineStarts = BooleanField("Do NOT add '--no_refine_starts' as option to TransDecoder.")
    outputFormat = SelectMultipleField(
        "Output Format"
        ,choices=[
            ("1","TSV Format")
            ,("2","CSV Format")
            ,("3","FASTA Amino Acid")
            ,("4","FASTA Nucleotide")
        ]
        ,description="Specify the output format for the processed alignments. Multiple formats can "
                     "be selected."
    )
    taxonomy = StringField(
        "Taxonomy"
        ,description="Specify the type of species/taxonomy you are analyzing and would like hits "
                     "closer in taxonomic relevance to be favored (based on NCBI Taxonomic "
                     "Database)."
    )
    queryCoverage = DecimalField(
        "Query Coverage"
        ,description="Specify the minimum query coverage to be allowed during similarity searching."
    )
    targetCoverage = DecimalField(
        "Target Coverage"
        ,description="Specify the minimum target coverage to be allowed during similarity "
                     "searching."
    )
    eValue = DecimalField(
        "E-Value"
        ,places=6
        ,description="Specify the E-Value that will be used as a cutoff during similarity "
                     "searching."
    )
    submit = SubmitField("Update")
