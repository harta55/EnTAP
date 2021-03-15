"""
Contains the GunzipTask class.
"""
from ..Abstract.AbstractTask import *
from ..Model.InputsModel import *
from flask import render_template
from os.path import join as pathJoin
from subprocess import run as pRun




class GunzipTask(AbstractTask):
    """
    This is the gunzip task class. It provides a task for running gunzip on a
    list of input files. It uses the path attribute from the inputs model to
    determine the full path of each file gunzipped, and therefore can only be
    used on input files.
    """


    def __init__(
        self
        ,fileNames
    ):
        """
        Initializes this new gunzip task with the given list of input files.

        Parameters
        ----------
        fileNames : list
                    Filenames of input files that this task will run gunzip on.
        """
        super().__init__()
        self.__fileNames = fileNames
        self._setRenderVars_(stage="init")


    def render(
        self
        ,**kwargs
    ):
        return render_template("task/gunzip.html",**kwargs)


    def run(
        self
    ):
        try:
            worked = []
            failed = []
            i = 1
            t = len(self.__fileNames)
            for fileName in self.__fileNames:
                self._setRenderVars_(stage="gunzip",fileName=fileName,current=i,total=t)
                path = pathJoin(InputsModel.PATH,fileName)
                cmd = ["gunzip",path]
                if pRun(cmd).returncode == 0:
                    worked.append(fileName)
                else:
                    failed.append(fileName)
                i += 1
            self._setRenderVars_(stage="fin",worked=worked,failed=failed)
            return not failed
        except Exception as e:
            self._setRenderVars_(stage="error",error=str(e))
            return False


    def title(
        self
    ):
        return "Gunzip"
