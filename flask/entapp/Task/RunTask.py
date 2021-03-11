"""
Contains the RunTask class.
"""
from ..Abstract.AbstractTask import *
from ..Controller.ConfigController import *
from ..Model.ConfigModel import *
from ..Model.DatabasesModel import *
from ..Model.LogsModel import *
from flask import render_template
from subprocess import run as pRun




class RunTask(AbstractTask):
    """
    This is the run task. It provides a task for running EnTAP from a given
    input FASTA file and optional BAM file if frame selection is enabled.
    """


    def __init__(
        self
        ,inputPath
        ,frameSelect=False
        ,bamPath=""
    ):
        """
        Initializes this new run task with the given input FASTA file, frame
        selection state, and BAM file if frame selection is enabled.

        Parameters
        ----------
        inputPath : string
                    Full absolute path of the input FASTA file.
        frameSelect : bool
                      True to enable frame selection or false otherwise.
        bamPath : string
                  Full absolute path of the input BAM file. Frame selection must
                  be enabled for this to be used.
        """
        super().__init__()
        self.__input = inputPath
        self.__fs = frameSelect
        self.__bam = bamPath
        self._setRenderVars_(stage="init")


    def render(
        self
        ,**kwargs
    ):
        if kwargs["stage"] == "running":
            output = ""
            logs = LogsModel()
            log = logs.newest()
            if log:
                output = log.tail()
            return render_template("task/setup.html",output=output,**kwargs)
        else:
            return render_template("task/setup.html",**kwargs)


    def run(
        self
    ):
        self._setRenderVars_(stage="running",frameSelect=self.__fs)
        config = ConfigModel()
        databases = DatabasesModel()
        cmd = ["EnTAP"]
        cmd.append("--runP" if self.__fs else "--runN")
        cmd += ["-i",self.__input]
        for dp in databases.enabledPaths():
            cmd += ["-d",dp]
        if self.__fs:
            cmd += ["-a",self.__bam]
        cmd += [
            "--ini"
            ,ConfigController.PATH
            ,"-t"
            ,str(config.threadNum())
            ,"--out-dir"
            ,DatabasesModel.OUT_PATH
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
        return "Run"
