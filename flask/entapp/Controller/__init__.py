"""
Contains controller classes and functions.
"""
from .ConfigController import *
from .TaskController import *
import flask
from flask import redirect
from flask import url_for


def update():
    """
    Updates the singleton task controller.

    Returns
    -------
    result : object
             A redirect to the status page if a task is currently running or
             else none.
    """
    if (
        not taskController.isSetup()
        and (
            flask.request.endpoint.startswith("DatabasesView:")
            or flask.request.endpoint.startswith("InputsView:")
            or flask.request.endpoint.startswith("RunView:")
        )
    ):
        return redirect(url_for("RootView:setup"))
    taskController.update()
    if taskController.isRunning() and flask.request.endpoint != "RootView:status":
        return redirect(url_for("RootView:status"))


taskController = TaskController()
configController = ConfigController()
