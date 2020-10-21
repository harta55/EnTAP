"""
Contains the BasicConfigForm class.
"""
from decimal import Decimal
import flask_wtf
import json
import wtforms as wtf
import wtforms.fields.html5 as wtf5








class BasicConfigForm(flask_wtf.FlaskForm):
    """
    Detailed description.
    """
    JSON_PATH = "/workspace/flask/basic_config.json"
    dataType = wtf.SelectField(
        "Data Type"
        ,choices=[("0","Serialized Database"),("1","SQLite Database")]
        ,description="Specifies which database EnTAP will download/generate."
    )
    fpkm = wtf5.DecimalField(
        "FPKM Threshold"
        ,description="Specifies the FPKM threshold with expression analysis. EnTAP will filter out "
                     "transcripts below this value!"
    )
    singleEnd = wtf.BooleanField("BAM/SAM files generated through single-end reads.")
    seqComplete = wtf.BooleanField("All sequences are complete proteins.")
    minProteinLength = wtf5.IntegerField(
        "Minimum Protein Length"
        ,description="Specifies the minimum protein length."
    )
    refineStarts = wtf.BooleanField("Do NOT add '--no_refine_starts' as option to TransDecoder.")
    outputFormat = wtf.SelectMultipleField(
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
    taxonomy = wtf.StringField(
        "Taxonomy"
        ,description="Specify the type of species/taxonomy you are analyzing and would like hits "
                     "closer in taxonomic relevance to be favored (based on NCBI Taxonomic "
                     "Database)."
    )
    queryCoverage = wtf5.DecimalField(
        "Query Coverage"
        ,description="Specify the minimum query coverage to be allowed during similarity searching."
    )
    targetCoverage = wtf5.DecimalField(
        "Target Coverage"
        ,description="Specify the minimum target coverage to be allowed during similarity "
                     "searching."
    )
    eValue = wtf5.DecimalField(
        "E-Value"
        ,places=6
        ,description="Specify the E-Value that will be used as a cutoff during similarity "
                     "searching."
    )
    submit = wtf.SubmitField("Update")


    def fromJson(
        self
        ,data
        ):
        """
        Detailed description.

        Parameters
        ----------
        data : object
               Detailed description.
        """
        data = json.loads(data)
        self.dataType.data = data["dataType"]
        self.fpkm.data = Decimal(data["fpkm"])
        self.singleEnd.data = data["singleEnd"]
        self.seqComplete.data = data["seqComplete"]
        self.minProteinLength.data = data["minProteinLength"]
        self.refineStarts.data = data["refineStarts"]
        self.outputFormat.data = data["outputFormat"]
        self.taxonomy.data = data["taxonomy"]
        self.queryCoverage.data = Decimal(data["queryCoverage"])
        self.targetCoverage.data = Decimal(data["targetCoverage"])
        self.eValue.data = Decimal(data["eValue"])


    def load(
        self
        ):
        """
        Detailed description.
        """
        with open(self.JSON_PATH,"r") as ifile:
            self.fromJson(ifile.read())


    def save(
        self
        ):
        """
        Detailed description.
        """
        with open(self.JSON_PATH,"w") as ofile:
            ofile.write(self.toJson()+"\n")


    def toJson(
        self
        ):
        """
        Detailed description.
        """
        data = {
            "dataType": self.dataType.data
            ,"fpkm": str(self.fpkm.data)
            ,"singleEnd": self.singleEnd.data
            ,"seqComplete": self.seqComplete.data
            ,"minProteinLength": self.minProteinLength.data
            ,"refineStarts": self.refineStarts.data
            ,"outputFormat": self.outputFormat.data
            ,"taxonomy": self.taxonomy.data
            ,"queryCoverage": str(self.queryCoverage.data)
            ,"targetCoverage": str(self.targetCoverage.data)
            ,"eValue": str(self.eValue.data)
        }
        return json.dumps(data,indent=4)
