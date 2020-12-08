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
        self.__url = url
        self.__total = 1
        self.__progress = 0
        self.__result = enums.TaskResult.Running


    def output(
        self
        ):
        """
        Detailed description.
        """
        return "Downloading "+self.__url+"\n"+self.__progress+"/"+self.__total+" Complete"


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
            rq = requests.get(url,stream=True)
            self.__total = rq.headers.get("content-length")
            for chunk in rq.iter_content(chunk_size=1024):
                if chunk:
                    self.__progress += 1024
            self.__result = enums.TaskResult.Finished
        except:
            self.__result = enums.TaskResult.Error
