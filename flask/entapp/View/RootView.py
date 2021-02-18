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


    def about(
        self
    ):
        """
        Detailed description.
        """
        return render_template("about.html")


    def index(
        self
    ):
        """
        Detailed description.
        """
        return render_template("index.html")


    def status(
        self
    ):
        """
        Detailed description.
        """
        return render_template("status.html",task=taskController)
