"""
Contains the ContamsConfigForm class.
"""
import json
import flask
import flask_wtf
import wtforms as wtf








class ContamsConfigForm(flask_wtf.FlaskForm):
    """
    Detailed description.
    """
    JSON_PATH = "/workspace/flask/contams_config.json"
    contams = wtf.SelectMultipleField(
        "Contaminants"
        ,choices=[]
        ,description="Select any number of contaminants to remove them."
    )
    rmSubmit = wtf.SubmitField("Remove Selected Contaminants")
    newContam = wtf.StringField(
        "New Contaminant"
        ,description="Use this field to input a new contaminant."
    )
    addSubmit = wtf.SubmitField("Add New Contaminant")


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
        self.contams.choices = [(c,c) for c in data]


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
        if self.submit:
            name = self.name.data
            if name in self.__contams:
                flash("New Contaminant already exists.","danger")
            else:
                self.__contams.add(name)
                flash("New Contaminant successfully added.","success")


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
        data = [c for (c,c) in self.contams.choices]
        return json.dumps(data,indent=4)
