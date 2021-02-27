"""
Contains the ContamsView class.
"""
from ..Controller import configController
from ..Form.ContamsForm import *
from ..Model.ContamsModel import *
from flask import flash
from flask import redirect
from flask import render_template
from flask import url_for
from flask_classful import FlaskView
from flask_classful import route




class ContamsView(FlaskView):
    """
    This is the contaminants view class. It provides a flask view for the
    contaminants list configuration of EnTAP.
    """
    route_base = "/contams/"


    @route("/add/",methods=["POST"])
    def add(
        self
    ):
        """
        Adds the new contaminant provided in the submitted form if the form
        passes all validation.

        Returns
        -------
        result : object
                 Flask redirect to this view's index page.
        """
        form = ContamsForm()
        if form.validate():
            contams = ContamsModel()
            name = form.name.data
            if name in contams:
                flash("Contaminant with given name already exists.","danger")
            else:
                contams.add(name)
                contams.save()
                configController.update()
                flash("New contaminant successfully added.","success")
        else:
            flash("Failed adding new contaminant.","danger")
        return redirect(url_for("ContamsView:index"))


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
        form = ContamsForm()
        contams = list(ContamsModel())
        contams.sort()
        return render_template("contams.html",contams=contams,form=form)


    def remove(
        self
        ,name
    ):
        """
        Removes the given contaminant.

        Parameters
        ----------
        name : string
               Name of removed contaminant.

        Returns
        -------
        result : object
                 Flask redirect to this view's index page.
        """
        contams = ContamsModel()
        if name in contams:
            contams.remove(name)
            contams.save()
            configController.update()
            flash("Contaminant successfully removed.","success")
        else:
            flash("Contaminant with given name does not exist.","danger")
        return redirect(url_for("ContamsView:index"))
