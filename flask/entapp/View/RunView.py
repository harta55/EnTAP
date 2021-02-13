"""
Contains the RunView class.
"""
from ..Controller import taskController
from ..Form.RunForm import *
from ..Model.DatabasesModel import *
from ..Model.InputsModel import *
from ..Task.RunTask import *
from flask import flash
from flask import redirect
from flask import render_template
from flask import url_for
from flask_classful import FlaskView
from flask_classful import route
from os.path import join as pathJoin




class RunView(FlaskView):
    """
    Detailed description.
    """
    route_base = "/run/"


    def index(
        self
    ):
        """
        Detailed description.
        """
        databases = DatabasesModel()
        inputs = InputsModel()
        form = RunForm([inputs.name(i) for i in inputs])
        return render_template("run.html",form=form,databases=databases.enabledNames())


    @route("/start/",methods=["POST"])
    def start(
        self
    ):
        """
        Detailed description.
        """
        inputs = InputsModel()
        form = RunForm([inputs.name(i) for i in inputs])
        if form.validate():
            taskController.start(RunTask(pathJoin(InputsModel.PATH,form.inputSelect.data)))
            flash("Started run task.","success")
            return redirect(url_for("RootView:status"))
        else:
            flash("Failed starting to run task.","danger")
            return redirect(url_for("RunView:index"))
