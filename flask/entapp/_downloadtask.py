"""
Contains the DownloadTask class.
"""
from . import enums
from . import interfaces
import os
import requests
from . import settings








class DownloadTask(interfaces.AbstractTask):
    """
    Detailed description.
    """


    def __init__(
        self
        ,url
        ):
        """
        Detailed description.

        Parameters
        ----------
        url : object
              Detailed description.
        """
        super().__init__()
        self.__url = url
        self.__fileName = url[url.rfind("/")+1:]
        self.__total = 1
        self.__progress = 0
        self.__result = enums.TaskResult.Running
        self.__finalOutput = None
        self.__error = ""


    def errorOutput(
        self
        ):
        """
        Detailed description.
        """
        return self.__error


    def finalOutput(
        self
        ):
        """
        Detailed description.
        """
        return self.__finalOutput


    def output(
        self
        ):
        """
        Detailed description.
        """
        ret1 = "Downloading "+self.__fileName+"\n"
        ret2 = None
        if not self.__total:
            ret2 = self.__reportBySize_()+" Downloaded"
        else:
            ret2 = ""
        return ret1+ret2


    def result(
        self
        ):
        """
        Detailed description.
        """
        return self.__result


    def run(
        self
        ):
        """
        Detailed description.
        """
        self.__progress = 0
        try:
            rq = requests.get(self.__url,stream=True)
            self.__total = int(rq.headers.get("Content-length",0))
            path = os.path.join(settings.DB_PATH,self.__fileName)
            if not path.endswith(".fa.gz"):
                self.__error = "Unknown file extension, remote file must end in '.fa.gz'."
                self.__result = enums.TaskResult.Error
                return
            elif os.path.exists(path) or os.path.exists(path[:-3]):
                self.__error = (
                    "'"
                    + self.__fileName[:-6]
                    + "' already exists as a database. Please remove it first if you wish to"
                      " overwrite it."
                )
                self.__result = enums.TaskResult.Error
                return
            with open(path,"wb") as ofile:
                for chunk in rq.iter_content(chunk_size=1024):
                    if chunk:
                        ofile.write(chunk)
                        self.__progress += 1024
            # GUNZIP...
            self.__result = enums.TaskResult.Finished
            self.__finalOutput = "Finished downloading "+self.__fileName
        except requests.exceptions.RequestException as e:
            self.__error = str(e)
            self.__result = enums.TaskResult.Error


    def title(
        self
        ):
        """
        Detailed description.
        """
        return "Remote Database Download"


    def __reportBySize_(
        self
        ):
        """
        Detailed description.
        """
        scale = ("B","KB","MB","GB","TB")
        size = self.__progress
        i = 0
        while size > 1024 and i < len(scale):
            size /= 1024
            i += 1
        return f"{size:.2f}"+scale[i]
