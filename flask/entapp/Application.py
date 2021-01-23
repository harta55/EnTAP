"""
Contains the Application class.
"""
from . import CONFIG_PATH
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
        self.__initConfig_()
        ContamsModel.initialize()


    def __initConfig_(
        self
    ):
        """
        Detailed description.
        """
        if not pathExists(dirname(CONFIG_PATH)):
            makedirs(dirname(CONFIG_PATH))
        if not pathExists(CONFIG_PATH):
            with open(CONFIG_PATH,"w") as ofile:
                ofile.write(
                    dumps(
                        {
                            "dataType": "1"
                            ,"fpkm": "0.50"
                            ,"singleEnd": True
                            ,"seqComplete": True
                            ,"minProteinLength": 100
                            ,"refineStarts": True
                            ,"outputFormat": [
                                "1"
                                ,"3"
                                ,"4"
                            ]
                            ,"taxonomy": ""
                            ,"queryCoverage": "50.00"
                            ,"targetCoverage": "50.00"
                            ,"eValue": "0.000001"
                        }
                            ,indent=4
                    )
                )
