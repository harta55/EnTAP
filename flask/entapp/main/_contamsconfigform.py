"""
Contains the ContamsConfigForm class.
"""
import flask_wtf
import wtforms as wtf








class ContamsConfigForm(flask_wtf.FlaskForm):
    """
    Detailed description.
    """
    contams = wtf.SelectMultipleField(
        "Contaminants"
        ,choices=[]
        ,description="Select any number of contaminants to remove them."
    )


    def fromJson(
        self
        ,data
        ):
        """
        Detailed description.

        Parameters
        ----------
        data : object
               Detailed description.
        """
        pass


    def load(
        self
        ):
        """
        Detailed description.
        """
        pass


    def save(
        self
        ):
        """
        Detailed description.
        """
        pass


    def toJson(
        self
        ):
        """
        Detailed description.
        """
        pass
