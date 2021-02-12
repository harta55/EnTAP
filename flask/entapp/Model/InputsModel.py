"""
Contains the InputsModel class.
"""
from os import listdir
from os import makedirs
from os import remove as rmFile
from os.path import exists as pathExists
from os.path import join as pathJoin




class InputsModel():
    """
    Detailed description.
    """
    PATH = "/workspace/entap/infiles"


    def __init__(
        self
    ):
        """
        Detailed description.
        """
        self.__files = listdir(self.PATH)


    def __contains__(
        self
        ,name
    ):
        """
        Detailed description.

        Parameters
        ----------
        name : 
        """
        return name in self.__files


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
        elif (
            ext == "fasta"
            or ext == "fna"
            or ext == "ffn"
            or ext == "faa"
            or ext == "fa"
            or ext == "frn"
        ):
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


    def remove(
        self
        ,names
    ):
        """
        Detailed description.

        Parameters
        ----------
        names : 
        """
        for name in names:
            if name not in self.__files:
                return False
            rmFile(pathJoin(self.PATH,name))
        return True
