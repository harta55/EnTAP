"""
Contains the UninformsForm class.
"""
from flask_wtf import FlaskForm
from wtforms import StringField
from wtforms import SubmitField




class UninformsForm(FlaskForm):
    """
    Detailed description.
    """
    name = StringField(
        "New Uninformative"
        ,description="Use this field to input a new uninformative."
    )
    submit = SubmitField("Add New Uninformative")
