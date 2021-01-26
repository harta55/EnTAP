"""
Contains the Application class.
"""
from .Model.ConfigModel import *
from .Model.ContamsModel import *
from flask import Flask
from json import dumps
from os.path import dirname
from os.path import exists as pathExists




class Application(Flask):
    """
    Detailed description.
    """


    def __init__(
        self
        ,*args
        ,**kwargs
    ):
        """
        Detailed description.

        Parameters
        ----------
        *args : 
        **kwargs : 
        """
        super().__init__(*args,**kwargs)
        ConfigModel.initialize()
        ContamsModel.initialize()
