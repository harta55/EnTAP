"""
Contains the IndexTask class.
"""
from ..Abstract.AbstractTask import *
from ..Controller.ConfigController import *
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
        super().__init__()
        self.__fileNames = fileNames
        self.__failed = []
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
            indexed = self.__index_(fileNames)
            self._setRenderVars_(
                stage="fin"
                ,indexed=indexed
                ,failed=self.__failed
            )
            return not self.__failed
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
            if pRun(cmd).returncode == 0:
                ret.append(fileName[:-3])
            else:
                self.__failed.append((fileName,"Failed running gunzip on file."))
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
        ret = []
        i = 0
        t = len(fileNames)
        for fileName in fileNames:
            path = pathJoin(DatabasesModel.PATH,fileName)
            cmd = [
                "EnTAP"
                ,"--config"
                ,"-d"
                ,path
                ,"--out-dir"
                ,DatabasesModel.OUT_PATH
                ,"-t"
                ,"8"
                ,"--ini"
                ,ConfigController.PATH
            ]
            self._setRenderVars_(stage="index",fileName=fileName,current=i,total=t)
            if pRun(cmd).returncode == 0:
                ret.append(fileName)
            else:
                self.__failed.append((fileName,"Failed indexing database."))
            i += 1
        return ret
