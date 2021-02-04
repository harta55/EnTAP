"""
Detailed description.
"""
from .. import application as app
from .TaskController import *
import flask
from flask import redirect
from flask import url_for


@app.before_request
def update():
    """
    Detailed description.
    """
    taskController.update()
    if taskController.isRunning() and flask.request.endpoint != "RootView:status":
        return redirect(url_for("RootView:status"))


taskController = TaskController()
