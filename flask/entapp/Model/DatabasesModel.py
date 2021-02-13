"""
Contains the DatabasesModel class.
"""
from json import dumps
from json import loads
from os import listdir
from os import makedirs
from os import remove as rmFile
from os.path import dirname
from os.path import exists as pathExists
from os.path import join as pathJoin




class DatabasesModel():
    """
    Detailed description.
    """
    PATH = "/workspace/flask/db"
    OUT_PATH = "/workspace/entap/outfiles"
    ENABLED_PATH = "/workspace/flask/databases.json"


    def __init__(
        self
    ):
        """
        Detailed description.
        """
        self.__files = listdir(self.PATH)
        with open(self.ENABLED_PATH,"r") as ifile:
            self.__enabled = loads(ifile.read())


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


    def disable(
        self
        ,name
    ):
        """
        Detailed description.

        Parameters
        ----------
        name : 
        """
        if name in self.__enabled:
            self.__enabled.remove(name)


    def enable(
        self
        ,name
    ):
        """
        Detailed description.

        Parameters
        ----------
        name : 
        """
        if name not in self.__enabled:
            self.__enabled.append(name)


    def enabledNames(
        self
    ):
        """
        Detailed description.
        """
        return self.__enabled


    def enabledPaths(
        self
    ):
        """
        Detailed description.
        """
        return [self.__indexPath_(f) for f in self.__enabled]


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
        if not pathExists(dirname(cls.ENABLED_PATH)):
            makedirs(dirname(cls.ENABLED_PATH))
        if not pathExists(cls.ENABLED_PATH):
            with open(cls.ENABLED_PATH,"w") as ofile:
                ofile.write(dumps([]))


    def isEnabled(
        self
        ,index
    ):
        """
        Detailed description.

        Parameters
        ----------
        index : 
        """
        return self.__files[index] in self.__enabled


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
        return pathExists(self.__indexPath_(self.__files[index]))


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
            indexPath = self.__indexPath_(name)
            if pathExists(indexPath):
                rmFile(indexPath)
        return True


    def save(
        self
    ):
        """
        Detailed description.
        """
        with open(self.ENABLED_PATH,"w") as ofile:
            ofile.write(dumps(self.__enabled))


    @classmethod
    def __indexPath_(
        cls
        ,name
    ):
        """
        Detailed description.

        Parameters
        ----------
        name : 
        """
        return pathJoin(cls.OUT_PATH,"bin",name[:name.rfind(".")]+".dmnd")
