"""
Contains the RootView class.
"""
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
        return render_template('index.html')
