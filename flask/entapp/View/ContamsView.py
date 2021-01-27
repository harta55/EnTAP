"""
Contains the ContamsView class.
"""
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
    Detailed description.
    """
    route_base = "/contams/"


    @route("/add/",methods=["POST"])
    def add(
        self
    ):
        """
        Detailed description.
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
                flash("New contaminant successfully added.","success")
        else:
            flash("Failed adding new contaminant.","danger")
        return redirect(url_for("ContamsView:index"))


    def index(
        self
    ):
        """
        Detailed description.
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
        Detailed description.

        Parameters
        ----------
        name : 
        """
        contams = ContamsModel()
        if name in contams:
            contams.remove(name)
            contams.save()
            flash("Contaminant successfully removed.","success")
        else:
            flash("Contaminant with given name does not exist.","danger")
        return redirect(url_for("ContamsView:index"))
