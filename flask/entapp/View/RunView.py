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
    This is the run view. It provides a flask view for running EnTAP itself.
    """
    route_base = "/run/"


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
        databases = DatabasesModel()
        inputs = InputsModel()
        form = RunForm([inputs.name(i) for i in inputs])
        return render_template("run.html",form=form,databases=databases.enabledNames())


    @route("/start/",methods=["POST"])
    def start(
        self
    ):
        """
        Starts a new run task of EnTAP with the submitted form's values if it
        passes validation.

        Returns
        -------
        result : object
                 Flask redirect to to the root view's status route if EnTAP was
                 successfully started otherwise this view's index page.
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
