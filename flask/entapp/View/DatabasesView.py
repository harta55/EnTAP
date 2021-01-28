"""
Contains the DatabasesView class.
"""
from ..Model.DatabasesModel import *
from flask import render_template
from flask_classful import FlaskView




class DatabasesView(FlaskView):
    """
    Detailed description.
    """
    route_base = "/databases/"


    def index(
        self
    ):
        """
        Detailed description.
        """
        databases = DatabasesModel()
        return render_template("databases/index.html",databases=databases)
