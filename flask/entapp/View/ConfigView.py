"""
Contains the ConfigView class.
"""
from ..Form.ConfigForm import *
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
        form.load()
        return render_template('config/basic.html',form=form)


    @route('/update/',methods=["POST"])
    def update(
        self
    ):
        """
        Detailed description.
        """
        form = ConfigForm()
        if form.validate():
            form.save()
            flash("Configuration successfully updated.","success")
            return redirect(url_for("ConfigView:index"))
        else:
            flash("Configuration failed to update.","danger")
            return redirect(url_for("ConfigView:index"))