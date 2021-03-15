"""
Contains the ConfigView class.
"""
from ..Controller import configController
from ..Form.ConfigForm import *
from ..Model.ConfigModel import *
from flask import flash
from flask import redirect
from flask import render_template
from flask import url_for
from flask_classful import FlaskView
from flask_classful import route




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
