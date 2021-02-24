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
    This is the inputs model class. It provides a model for all input FASTA
    files contained in this flask application's local environment.
    """
    PATH = "/workspace/entap/infiles"


    def __init__(
        self
    ):
        self.__files = listdir(self.PATH)


    def __contains__(
        self
        ,name
    ):
        return name in self.__files


    def __len__(
        self
    ):
        return len(self.__files)


    def __iter__(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : range_iterator
                 A ranged iterator that will iterate through every valid integer
                 index for all input files of this model, starting at 0 and
                 ending at the last input file in the list.
        """
        return range(len(self.__files)).__iter__()


    def extName(
        self
        ,index
    ):
        """
        Getter method.

        Parameters
        ----------
        index : int
                Index of the input file whose extension name is returned.

        Returns
        -------
        result : string
                 Name of the extension of this model's input file at the given
                 index.
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
        Creates the directory where input files are stored if it does not exist.
        """
        if not pathExists(cls.PATH):
            makedirs(cls.PATH)


    def name(
        self
        ,index
    ):
        """
        Getter method.

        Parameters
        ----------
        index : int
                Index of the input file whose filename is returned.

        Returns
        -------
        result : string
                 Filename of this model's input file at the given index.
        """
        return self.__files[index]


    def remove(
        self
        ,names
    ):
        """
        Removes the given list of input files from this model.

        Parameters
        ----------
        names : list
                Filenames of input files that are removed from this model.
        """
        for name in names:
            if name not in self.__files:
                return False
            rmFile(pathJoin(self.PATH,name))
        return True
