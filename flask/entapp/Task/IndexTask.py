"""
Contains the IndexTask class.
"""
from ..Abstract.AbstractTask import *
from ..Model.DatabasesModel import *
from flask import render_template
from os.path import join as pathJoin
from subprocess import run as pRun




class IndexTask(AbstractTask):
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
        self.__fileNames = fileNames
        self.__gErrors = []
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
        return render_template("task/index.html",**kwargs)


    def run(
        self
    ):
        """
        Detailed description.
        """
        try:
            fileNames = []
            gFileNames = []
            for fileName in self.__fileNames:
                if fileName.endswith(".gz"):
                    gFileNames.append(fileName)
                else:
                    fileNames.append(fileName)
            fileNames += self.__gunzip_(gFileNames)
            indexSuccess = self.__index_(fileNames)
            self._setRenderVars_(
                stage="fin"
                ,indexSuccess=indexSuccess
                ,indexedNames=fileNames
                ,gunzipErrors=self.__gErrors
            )
            return indexSuccess and not self.__gErrors
        except Exception as e:
            self._setRenderVars_(stage="error",error=str(e))
            return False


    def title(
        self
    ):
        """
        Detailed description.
        """
        return "Index"


    def __gunzip_(
        self
        ,fileNames
    ):
        """
        Detailed description.

        Parameters
        ----------
        fileNames : 
        """
        ret = []
        i = 0
        t = len(fileNames)
        for fileName in fileNames:
            path = pathJoin(DatabasesModel.PATH,fileName)
            cmd = ["gunzip",path]
            self._setRenderVars_(stage="gunzip",fileName=fileName,current=i,total=t)
            if pRun(cmd).returncode != 0:
                ret.append(fileName[:-3])
            else:
                self.__gErrors = (fileName,"Failed running gunzip on file.")
            i += 1
        return ret


    def __index_(
        self
        ,fileNames
    ):
        """
        Detailed description.

        Parameters
        ----------
        fileNames : 
        """
        return True
