"""
Contains the SetupTask class.
"""
from ..Abstract.AbstractTask import *
from ..Controller.ConfigController import *
from ..Model.ConfigModel import *
from ..Model.DatabasesModel import *
from ..Model.LogsModel import *
from flask import render_template
from subprocess import run as pRun




class SetupTask(AbstractTask):
    """
    This is the setup task. It provides a task for running EnTAP configuration
    for the first time, allowing it to download its root database and index it.
    """


    def __init__(
        self
    ):
        super().__init__()
        self._setRenderVars_(stage="init")


    def render(
        self
        ,**kwargs
    ):
        return render_template("task/setup.html",**kwargs)


    def run(
        self
    ):
        logs = LogsModel()
        self._setRenderVars_(stage="running",output=logs.newest().tail())
        config = ConfigModel()
        cmd = [
            "EnTAP"
            ,"--config"
            ,"--out-dir"
            ,DatabasesModel.OUT_PATH
            ,"-t"
            ,str(config.threadNum())
            ,"--ini"
            ,ConfigController.PATH
        ]
        if pRun(cmd).returncode == 0:
            self._setRenderVars_(stage="success")
            return True
        else:
            self._setRenderVars_(stage="failed")
            return False


    def title(
        self
    ):
        return "Setup"
