"""
Contains the Application class.
"""
from .Controller import configController
from .Model.ConfigModel import *
from .Model.ContamsModel import *
from .Model.DatabasesModel import *
from .Model.UninformsModel import *
from flask import Flask
from json import dumps




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
        DatabasesModel.initialize()
        UninformsModel.initialize()
        configController.update()
