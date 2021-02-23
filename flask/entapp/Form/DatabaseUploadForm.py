"""
Contains the DatabaseUploadForm class.
"""
from flask_wtf import FlaskForm
from wtforms import TextAreaField
from wtforms import SelectMultipleField
from wtforms import SubmitField




class DatabaseUploadForm(FlaskForm):
    """
    This is the database upload form. It provides input for remotely uploading
    databases from predetermined URLs and custom URLs.
    """
    defUrls = SelectMultipleField(
        "Databases"
        ,choices=[
            (
                "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
                ,"NCBI Non-Redundant (NR)"
            )
            ,(
                "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/"
                "complete/uniprot_sprot.fasta.gz"
                ,"ExPASy Swiss-Prot"
            )
            ,(
                "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/"
                 "complete/uniprot_trembl.fasta.gz"
                ,"ExPASy TrEMBL"
            )
        ]
        ,description="Select any number of databases that EnTAP will download to its database"
                     " folder one by one as a task."
    )
    urls = TextAreaField(
        "Remote URLs"
        ,description="Specify a list of remote URLs, separated by one URL per line, that EnTAP will"
                     " download to its database folder one by one as a task."
    )
    submit = SubmitField("Start Remote Database Upload Task")
    cancel = SubmitField("Cancel")
