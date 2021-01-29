"""
Contains the DatabasesModel class.
"""
from os import listdir
from os import makedirs
from os.path import exists as pathExists




class DatabasesModel():
    """
    Detailed description.
    """
    PATH = "/workspace/flask/db"


    def __init__(
        self
    ):
        """
        Detailed description.
        """
        self.__files = listdir(self.PATH)


    def __len__(
        self
    ):
        """
        Detailed description.
        """
        return len(self.__files)


    def __iter__(
        self
    ):
        """
        Detailed description.
        """
        return range(len(self.__files)).__iter__()


    def extName(
        self
        ,index
    ):
        """
        Detailed description.

        Parameters
        ----------
        index : 
        """
        name = self.__files[index]
        ext = name[name.rfind(".")+1:]
        if ext == "gz":
            return "GUNZIP"
        elif ext == "fa" or ext == "faa" or ext == "fasta":
            return "FASTA"
        else:
            return "UNKNOWN"


    @classmethod
    def initialize(
        cls
    ):
        """
        Detailed description.
        """
        if not pathExists(cls.PATH):
            makedirs(cls.PATH)


    def isIndexed(
        self
        ,index
    ):
        """
        Detailed description.

        Parameters
        ----------
        index : 
        """
        # TODO
        return False


    def name(
        self
        ,index
    ):
        """
        Detailed description.

        Parameters
        ----------
        index : 
        """
        return self.__files[index]
