"""
Detailed description.
"""
from .Application import *
from .View.ConfigView import *
from .View.ContamsView import *
from .View.RootView import *
from flask import Flask


def createApp():
    """
    Detailed description.
    """
    app = Application(__name__)
    app.config["SECRET_KEY"] = "hard to guess phrase DO NOT USE IN PRODUCTION"
    ConfigView.register(app)
    ContamsView.register(app)
    RootView.register(app)
    return app
