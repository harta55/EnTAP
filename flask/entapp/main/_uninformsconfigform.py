"""
Contains the UninformsConfigForm class.
"""
import json
import flask
import flask_wtf
import wtforms as wtf








class UninformsConfigForm(flask_wtf.FlaskForm):
    """
    Detailed description.
    """
    JSON_PATH = "/workspace/flask/uninforms_config.json"
    uninforms = wtf.SelectMultipleField(
        "Uninformatives"
        ,choices=[]
        ,description="Select any number of uninformatives to remove them."
    )
    rmSubmit = wtf.SubmitField("Remove Selected Uninformatives")
    newUninform = wtf.StringField(
        "New Uninformative"
        ,description="Use this field to input a new uninformative."
    )
    addSubmit = wtf.SubmitField("Add New Uninformative")


    def __init__(
        self
        ):
        """
        Detailed description.
        """
        super().__init__()
        self.load()


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
        data = json.loads(data)
        self.uninforms.choices = [(c,c) for c in data]


    def load(
        self
        ):
        """
        Detailed description.
        """
        with open(self.JSON_PATH,"r") as ifile:
            self.fromJson(ifile.read())


    def processContams(
        self
        ):
        """
        Detailed description.
        """
        if self.addSubmit and self.newUninform.data:
            c = self.newUninform.data
            c = (c,c)
            if not c in self.uninforms.choices:
                self.uninforms.choices.append(c)
                flask.flash("New Uninformative successfully added.","success")
            else:
                flask.flash("New Uninformative already exists.","danger")
        elif self.rmSubmit and self.uninforms.data:
            for r in self.uninforms.data:
                self.uninforms.choices.remove((r,r))
            self.uninforms.data.clear()
            flask.flash("Selected Uninformatives successfully removed.","success")


    def save(
        self
        ):
        """
        Detailed description.
        """
        with open(self.JSON_PATH,"w") as ofile:
            ofile.write(self.toJson()+"\n")


    def toJson(
        self
        ):
        """
        Detailed description.
        """
        data = [c for (c,c) in self.uninforms.choices]
        return json.dumps(data,indent=4)
