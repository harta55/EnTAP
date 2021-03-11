"""
Contains the RunForm class.
"""
from flask_wtf import FlaskForm
from wtforms import BooleanField
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
    frameSelection = BooleanField("Enable frame selection (only applicable to nucleotide sequences).")
    bamSelect = SelectField(
        "Select SAM/BAM"
        ,choices=[]
        ,description="Select a SAM/BAM file for processing from the list provided."
    )
    run = SubmitField("Start")


    def __init__(
        self
        ,inputs
        ,bams
    ):
        """
        Initializes this new run form with the given list of input and bam
        files.

        Parameters
        ----------
        inputs : list
                 Filenames of all valid input FASTA files.
        bams : list
               Filenames of all valid input BAM files.
        """
        super().__init__()
        for i in inputs:
            self.inputSelect.choices.append((i,i))
        self.bamSelect.choices.append(("","None"))
        for bam in bams:
            self.bamSelect.choices.append((bam,bam))
