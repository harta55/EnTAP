"""
Contains the RunTask class.
"""
from ..Abstract.AbstractTask import *
from ..Controller.ConfigController import *
from ..Model.ConfigModel import *
from ..Model.DatabasesModel import *
from flask import render_template
from subprocess import run as pRun




class RunTask(AbstractTask):
    """
    Detailed description.
    """


    def __init__(
        self
        ,path
    ):
        """
        Detailed description.

        Parameters
        ----------
        path : 
        """
        super().__init__()
        self.__path = path
        self._setRenderVars_(stage="init")


    def render(
        self
        ,**kwargs
    ):
        """
        Detailed description.

        Parameters
        ----------
        **kwargs : 
        """
        return render_template("task/run.html",**kwargs)


    def run(
        self
    ):
        """
        Detailed description.
        """
        config = ConfigModel()
        databases = DatabasesModel()
        cmd = ["EnTAP","--runP","-i",self.__path]
        for dp in databases.enabledPaths():
            cmd += ["-d",dp]
        cmd += [
            "--ini"
            ,ConfigController.PATH
            ,"-t"
            ,str(config.threadNum())
            ,"--out-dir"
            ,DatabasesModel.OUT_PATH
        ]
        print(cmd)
        if pRun(cmd).returncode == 0:
            self._setRenderVars_(stage="success")
            return True
        else:
            self._setRenderVars_(stage="failed")
            return False


    def title(
        self
    ):
        """
        Detailed description.
        """
        return "Run"
