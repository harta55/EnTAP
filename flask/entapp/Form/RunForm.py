"""
Contains the RunForm class.
"""
from flask_wtf import FlaskForm
from wtforms import SelectField
from wtforms import SubmitField




class RunForm(FlaskForm):
    """
    This is the run form. It provides all necessary and optional options for
    starting a EnTAP run job.
    """
    inputSelect = SelectField(
        "Select Input"
        ,choices=[]
        ,description="Select an input for processing from the list of inputs provided."
    )
    run = SubmitField("Start")


    def __init__(
        self
        ,names
    ):
        """
        Detailed description.

        Parameters
        ----------
        names : 
        """
        super().__init__()
        for name in names:
            self.inputSelect.choices.append((name,name))
