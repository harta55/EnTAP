"""
Contains the Application class.
"""
from .Controller import configController
from .Model.ConfigModel import *
from .Model.ContamsModel import *
from .Model.DatabasesModel import *
from .Model.InputsModel import *
from .Model.LogsModel import *
from .Model.UninformsModel import *
from flask import Flask
from json import dumps




class Application(Flask):
    """
    This is the EnTAP flask application class. It extends the flask application
    class, providing some custom initialization code this this specific
    application.
    """


    def __init__(
        self
        ,*args
        ,**kwargs
    ):
        """
        Initializes this new EnTAP flask application, padding all arguments onto
        the flask application's initialization.

        Parameters
        ----------
        *args : list
                All positional arguments.
        **kwargs : dictionary
                   All keyword arguments.
        """
        super().__init__(*args,**kwargs)
        ConfigModel.initialize()
        ContamsModel.initialize()
        DatabasesModel.initialize()
        InputsModel.initialize()
        LogsModel.initialize()
        UninformsModel.initialize()
        configController.update()
