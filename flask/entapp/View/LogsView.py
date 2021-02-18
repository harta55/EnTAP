"""
Contains the LogsView class.
"""
from ..Model.LogsModel import *
from flask import render_template
from flask_classful import FlaskView




class LogsView(FlaskView):
    """
    Detailed description.
    """
    route_base = "/logs/"


    def debug(
        self
        ,index
    ):
        """
        Detailed description.

        Parameters
        ----------
        index : 
        """
        index = int(index)
        logs = LogsModel()
        return render_template("logs/debug.html",log=logs[index])


    def index(
        self
    ):
        """
        Detailed description.
        """
        logs = LogsModel()
        return render_template("logs/index.html",logs=logs)


    def read(
        self
        ,index
    ):
        """
        Detailed description.

        Parameters
        ----------
        index : 
        """
        index = int(index)
        logs = LogsModel()
        return render_template("logs/read.html",log=logs[index])
