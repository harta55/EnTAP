"""
Contains the UninformsView class.
"""
from ..Controller import configController
from ..Form.UninformsForm import *
from ..Model.UninformsModel import *
from flask import flash
from flask import redirect
from flask import render_template
from flask import url_for
from flask_classful import FlaskView
from flask_classful import route




class UninformsView(FlaskView):
    """
    This is the contaminants view class. It provides a flask view for the
    uninformative list configuration of EnTAP.
    """
    route_base = "/uninforms/"


    @route("/add/",methods=["POST"])
    def add(
        self
    ):
        """
        Adds the new uninformative provided in the submitted form if the form
        passes all validation.

        Returns
        -------
        result : object
                 Flask redirect to this view's index page.
        """
        form = UninformsForm()
        if form.validate():
            uninforms = UninformsModel()
            name = form.name.data
            if not name:
                flash("Uninformative cannot be blank string.","danger")
            elif name in uninforms:
                flash("Uninformative with given name already exists.","danger")
            else:
                uninforms.add(name)
                uninforms.save()
                configController.update()
                flash("New uninformative successfully added.","success")
        else:
            flash("Failed adding new uninformative.","danger")
        return redirect(url_for("UninformsView:index"))


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
        form = UninformsForm()
        uninforms = list(UninformsModel())
        uninforms.sort()
        return render_template("uninforms.html",uninforms=uninforms,form=form)


    def remove(
        self
        ,name
    ):
        """
        Removes the given uninformative.

        Parameters
        ----------
        name : string
               Name of removed uninformative.

        Returns
        -------
        result : object
                 Flask redirect to this view's index page.
        """
        uninforms = UninformsModel()
        if name in uninforms:
            uninforms.remove(name)
            uninforms.save()
            configController.update()
            flash("Uninformative successfully removed.","success")
        else:
            flash("Uninformative with given name does not exist.","danger")
        return redirect(url_for("UninformsView:index"))
