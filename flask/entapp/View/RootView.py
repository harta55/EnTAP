"""
Contains the RootView class.
"""
from ..Controller import taskController
from flask import render_template
from flask_classful import FlaskView




class RootView(FlaskView):
    """
    Detailed description.
    """
    route_base = "/"


    def index(
        self
    ):
        """
        Detailed description.
        """
        return render_template("index.html")


    def run(
        self
    ):
        """
        Detailed description.
        """
        return render_template("run.html")


    def status(
        self
    ):
        """
        Detailed description.
        """
        return render_template("status.html",task=taskController)
