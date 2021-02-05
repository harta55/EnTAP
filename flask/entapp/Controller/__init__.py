"""
Detailed description.
"""
from .ConfigController import *
from .TaskController import *
import flask
from flask import redirect
from flask import url_for


def update():
    """
    Detailed description.
    """
    taskController.update()
    if taskController.isRunning() and flask.request.endpoint != "RootView:status":
        return redirect(url_for("RootView:status"))


taskController = TaskController()
configController = ConfigController()
