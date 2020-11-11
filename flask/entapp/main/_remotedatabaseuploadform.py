"""
Contains the RemoteDatabaseUploadForm class.
"""
import flask_wtf
import wtforms as wtf








class RemoteDatabaseUploadForm(flask_wtf.FlaskForm):
    """
    Detailed description.
    """
    customURL = wtf.StringField(
        "Custom URL"
        ,description="Specify a custom URL to download a gunzipped fasta database file."
    )
    addSubmit = wtf.SubmitField("Run Remote Upload")
