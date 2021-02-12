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
    Detailed description.
    """


    def __init__(
        self
        ,fileNames
    ):
        """
        Detailed description.

        Parameters
        ----------
        fileNames : 
        """
        super().__init__()
        self.__fileNames = fileNames
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
        return render_template("task/gunzip.html",**kwargs)


    def run(
        self
    ):
        """
        Detailed description.
        """
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
        """
        Detailed description.
        """
        return "Gunzip"
