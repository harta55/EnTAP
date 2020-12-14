"""
Contains the DownloadTask class.
"""
from . import enums
import requests
import threading








class DownloadTask(threading.Thread):
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


    def output(
        self
        ):
        """
        Detailed description.
        """
        if self.__finalOutput:
            return ("Remote Database Download",self.__finalOutput)
        else:
            ret1 = "Downloading "+self.__fileName+"\n"
            ret2 = None
            if not self.__total:
                ret2 = self.__reportBySize_()+" Downloaded"
            else:
                ret2 = ""
            return ("Remote Database Download",ret1+ret2)


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
            print(rq.headers)
            self.__total = int(rq.headers.get("Content-length",0))
            for chunk in rq.iter_content(chunk_size=1024):
                if chunk:
                    self.__progress += 1024
            self.__result = enums.TaskResult.Finished
            self.__finalOutput = "Finished downloading "+self.__finalName
        except:
            self.__result = enums.TaskResult.Error


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
