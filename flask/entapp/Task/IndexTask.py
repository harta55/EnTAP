"""
Contains the IndexTask class.
"""
from ..Abstract.AbstractTask import *
from ..Controller.ConfigController import *
from ..Model.ConfigModel import *
from ..Model.DatabasesModel import *
from flask import render_template
from os.path import join as pathJoin
from subprocess import run as pRun




class IndexTask(AbstractTask):
    """
    This is the index task class. It provides a task for running an EnTAP config
    command to index a list of databases, gunzipping any that have a gunzip
    extension. It uses the path attribute from the databases model to determine
    the full path of each database.
    """


    def __init__(
        self
        ,fileNames
    ):
        """
        Initializes this new index task with the given list of databases.

        Parameters
        ----------
        fileNames : list
                    Filenames of databases that this task will index, and gunzip
                    if necessary.
        """
        super().__init__()
        self.__fileNames = fileNames
        self.__failed = []
        self._setRenderVars_(stage="init")


    def render(
        self
        ,**kwargs
    ):
        return render_template("task/index.html",**kwargs)


    def run(
        self
    ):
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
        return "Index"


    def __gunzip_(
        self
        ,fileNames
    ):
        """
        Runs gunzip on the list of given databases.

        Parameters
        ----------
        fileNames : list
                    Filenames of databases that are gunzipped.

        Returns
        -------
        result : list
                 Filenames of databases that were successfully gunzipped. The
                 filenames are the new names without the gunzip extension.
        """
        ret = []
        i = 1
        t = len(fileNames)
        for fileName in fileNames:
            path = pathJoin(DatabasesModel.PATH,fileName)
            if pathExists(path[:-3]):
                self.__failed.append(
                    (fileName,"Unzipped version of file already exists, please remove it first.")
                )
                continue
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
        Runs EnTAP configuration on given list of databases one at a time.

        Parameters
        ----------
        fileNames : list
                    Filenames of databases that are indexed.

        Returns
        -------
        result : list
                 Filenames of databases that were successfully indexed by EnTAP.
        """
        ret = []
        i = 1
        t = len(fileNames)
        config = ConfigModel()
        databases = DatabasesModel()
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
                ,str(config.threadNum())
                ,"--ini"
                ,ConfigController.PATH
            ]
            self._setRenderVars_(stage="index",fileName=fileName,current=i,total=t)
            if pRun(cmd).returncode == 0:
                databases.enable(fileName)
                databases.save()
                ret.append(fileName)
            else:
                self.__failed.append((fileName,"Failed indexing database."))
            i += 1
        return ret
