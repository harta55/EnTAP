"""
Contains the ConfigView class.
"""
from ..Controller import configController
from ..Controller.ConfigController import *
from ..Model.DatabasesModel import *
from ..Form.ConfigForm import *
from ..Model.ConfigModel import *
from flask import flash
from flask import redirect
from flask import render_template
from flask import url_for
from flask_classful import FlaskView
from flask_classful import route
from json import dumps
from json import loads
from subprocess import check_output as cRun




class ConfigView(FlaskView):
    """
    This is the configuration view class. It provides a flask view for the basic
    configuration of EnTAP.
    """
    route_base = "/config/"


    def index(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : object
                 Index page of this view.
        """
        form = ConfigForm()
        config = ConfigModel()
        config.populate(form)
        return render_template("config.html",form=form)


    @route("/taxonomy/<tax>")
    def isTaxValid(
        self
        ,tax
    ):
        """
        Getter method.

        Parameters
        ----------
        tax : string
              A taxonomy field that is checked for its validity.

        Returns
        -------
        result : object
                 JSON encoded RESTful reply informing if the given taxonomy
                 string is valid. The JSON contains an object with the key
                 "isvalid" with a boolean value.
        """
        data = loads(
            cRun(
                [
                    "EnTAP"
                    ,"--api-taxon"
                    ,tax
                    ,"--out-dir"
                    ,DatabasesModel.OUT_PATH
                    ,"--ini"
                    ,ConfigController.PATH
                ]
            )
        )
        isValid = True if (data["error"] == False and data["valid"] == True) else False
        return dumps({"valid": isValid})


    @route("/update/",methods=["POST"])
    def update(
        self
    ):
        """
        Updates basic configuration with the values from the submitted form if
        they pass form validation.

        Returns
        -------
        result : object
                 Flask redirect back to this view's index page.
        """
        form = ConfigForm()
        if form.validate():
            config = ConfigModel()
            config.update(form)
            config.save()
            configController.update()
            flash("Configuration successfully updated.","success")
        else:
            flash("Configuration failed to update.","danger")
        return redirect(url_for("ConfigView:index"))
