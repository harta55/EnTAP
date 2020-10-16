"""
Detailed description.
"""
import flask_wtf
import wtforms as wtf
from wtforms import validators as fvs








class BasicConfigForm(flask_wtf.FlaskForm):
    """
    Detailed description.
    """
    dataType = wtf.SelectField(
        "Data Type"
        ,choices=[(0,"Serialized Database"),(1,"SQLite Database")]
        ,description="Specifies which database EnTAP will download/generate."
    )
    fpkm = wtf.DecimalField(
        "FPKM Threshold"
        ,description="Specifies the FPKM threshold with expression analysis. EnTAP will filter out "
                     "transcripts below this value!"
    )
    singleEnd = wtf.SelectField(
        "BAM/SAM File Generation"
        ,choices=[(False,"Paired End"),(True,"Single End")]
        ,description="Specifies if your BAM/SAM file was generated through single-end or "
                     "paired-end reads. This is only required in expression analysis."
    )
    submit = wtf.SubmitField("Update")
