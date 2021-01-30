"""
Contains the DatabaseUploadForm class.
"""
from flask_wtf import FlaskForm
from wtforms import StringField
from wtforms import SubmitField




class DatabaseUploadForm(FlaskForm):
    """
    Detailed description.
    """
    url = StringField(
        "Remote URL"
        ,description="Specify a remote URL that EnTAP will download to its database folder as a task."
    )
    submit = SubmitField("Start Remote Database Upload Task")
    cancel = SubmitField("Cancel")
