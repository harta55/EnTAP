"""
Contains the DatabaseUploadForm class.
"""
from flask_wtf import FlaskForm
from wtforms import TextAreaField
from wtforms import SubmitField




class DatabaseUploadForm(FlaskForm):
    """
    Detailed description.
    """
    url = TextAreaField(
        "Remote URLs"
        ,description="Specify a list of remote URLs, separated by one URL per line, that EnTAP will download to its database folder one by one as a task."
    )
    submit = SubmitField("Start Remote Database Upload Task")
    cancel = SubmitField("Cancel")
