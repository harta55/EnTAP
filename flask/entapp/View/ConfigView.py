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
    Detailed description.
    """
    route_base = "/config/"


    def index(
        self
    ):
        """
        Detailed description.
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
        Detailed description.
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
