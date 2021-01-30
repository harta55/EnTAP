"""
Contains the RemoteUploadTask class.
"""
from ..Abstract.AbstractTask import *
from ..Model.DatabasesModel import *
from requests import get as reqGet
from os.path import exists as pathExists
from os.path import join as pathJoin




class RemoteUploadTask(AbstractTask):
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
        url : 
        """
        super().__init__()
        self.__url = url


    def run(
        self
    ):
        """
        Detailed description.
        """
        try:
            fileName = self.__url[self.__url.rfind("/")+1:]
            progress = 0
            rq = reqGet(self.__url,stream=True)
            total = int(rq.headers.get("Content-length",0))
            path = pathJoin(DatabasesModel.PATH,fileName)
            if pathExists(path):
                self._clear_()
                self._print_(
                    "'"
                    + fileName
                    + "' already exists. Please remove it first if you wish to overwrite it."
                )
                self._flush_()
                return False
            with open(path,"wb") as ofile:
                for chunk in rq.iter_content(chunk_size=1024):
                    if chunk:
                        ofile.write(chunk)
                        progress += 1024
                        self._clear_()
                        self._print_("Downloading "+fileName+"\n")
                        if not total:
                            self._print_(self.__reportBySize_(progress)+" Downloaded")
                        else:
                            self._print_("%i% Complete"%(int(progress*100/total),))
                        self._flush_()
            self._clear_()
            self._print_("Finished downloading "+fileName)
            self._flush_()
            return True
        except Exception as e:
            self._clear_()
            self._print_(str(e))
            self._flush_()
            return False


    def title(
        self
    ):
        """
        Detailed description.
        """
        return "Remote Upload"


    @staticmethod
    def __reportBySize_(
        size
    ):
        """
        Detailed description.

        Parameters
        ----------
        size : 
        """
        scale = ("B","KB","MB","GB","TB")
        i = 0
        while size > 1024 and i < len(scale):
            size /= 1024
            i += 1
        return f"{size:.2f}"+scale[i]
