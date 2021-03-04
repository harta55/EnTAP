"""
Contains the RootView class.
"""
from ..Controller import taskController
from ..Task.SetupTask import *
from flask import flash
from flask import redirect
from flask import render_template
from flask import url_for
from flask_classful import FlaskView
from flask_classful import route




class RootView(FlaskView):
    """
    This is the root view class. It provides a flask view for basic the root
    pages of home, status, and about.
    """
    route_base = "/"


    def about(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : object
                 This view's about page.
        """
        return render_template("about.html")


    def index(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : object
                 This view's index home page.
        """
        return render_template("index.html")


    def setup(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : object
                 This view's setup page.
        """
        return render_template("setup.html")


    @route("/setup/start/")
    def startSetup(
        self
    ):
        """
        Runs the setup task for EnTAP.

        Returns
        -------
        result : object
                 Flask redirect to this view's status page.
        """
        taskController.start(SetupTask())
        flash("Started setup task.","success")
        return redirect(url_for("RootView:status"))


    def status(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : object
                 This view's status page.
        """
        return render_template("status.html",task=taskController)
