"""
Contains the ContamsForm class.
"""
from flask_wtf import FlaskForm
from wtforms import StringField
from wtforms import SubmitField




class ContamsForm(FlaskForm):
    """
    Detailed description.
    """
    name = StringField(
        "New Contaminant"
        ,description="Use this field to input a new contaminant."
    )
    submit = SubmitField("Add New Contaminant")
