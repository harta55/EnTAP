"""
Contains the GunzipTask class.
"""
from ..Abstract.AbstractTask import *
from ..Model.DatabasesModel import *
from os.path import join as pathJoin
from subprocess import run as pRun




class GunzipTask(AbstractTask):
    """
    Detailed description.
    """


    def __init__(
        self
        ,fileName
    ):
        """
        Detailed description.

        Parameters
        ----------
        fileName : 
        """
        super().__init__()
        self.__fileName = fileName


    def run(
        self
    ):
        """
        Detailed description.
        """
        try:
            path = pathJoin(DatabasesModel.PATH,self.__fileName)
            cmd = ["gunzip",path]
            self._clear_()
            self._print_("Gunzipping database '"+self.__fileName+"'.")
            self._flush_()
            assert(pRun(cmd).returncode==0)
            self._clear_()
            self._print_("Gunzipped database '"+self.__fileName+"'.")
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
        return "Gunzip"
