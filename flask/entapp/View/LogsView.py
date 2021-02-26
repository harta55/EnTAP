"""
Contains the LogsView class.
"""
from ..Model.LogsModel import *
from flask import render_template
from flask_classful import FlaskView




class LogsView(FlaskView):
    """
    This is the logs view. It provides a flask view for listing all EnTAP logs
    and viewing individual logs.
    """
    route_base = "/logs/"


    def debug(
        self
        ,index
    ):
        """
        Getter method.

        Parameters
        ----------
        index : int
                The index of the log whose full debug log is displayed.

        Returns
        -------
        result : object
                 This view's debug page, displaying the full debug log of the
                 EnTAP run at the given index.
        """
        index = int(index)
        logs = LogsModel()
        return render_template("logs/debug.html",log=logs[index])


    def index(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : object
                 This view's index page, listing all EnTAP logs.
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
        index : int
                The index of the log whose full log is displayed.

        Returns
        -------
        result : object
                 This view's read page, displaying the full log of the EnTAP run
                 at the given index.
        """
        index = int(index)
        logs = LogsModel()
        return render_template("logs/read.html",log=logs[index])
