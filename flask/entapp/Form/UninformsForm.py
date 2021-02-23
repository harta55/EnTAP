"""
Contains the UninformsForm class.
"""
from flask_wtf import FlaskForm
from wtforms import StringField
from wtforms import SubmitField




class UninformsForm(FlaskForm):
    """
    This is the uninformative configuration form. It provides input for adding a
    new uninformative to its respective list.
    """
    name = StringField(
        "New Uninformative"
        ,description="Use this field to input a new uninformative."
    )
    submit = SubmitField("Add New Uninformative")
