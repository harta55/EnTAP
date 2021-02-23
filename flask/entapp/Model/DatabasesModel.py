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
    This is the databases model class. It provides a model for all FASTA
    databases contained in this flask application's local environment. It also
    keeps track of important properties for each database. These properties are
    file extension, if it is indexed by EnTAP, and if it is enabled for use when
    running EnTAP.

    This database keeps track of enabled/disabled databases through its own
    custom JSON data file. Because of this the model must be saved when changes
    are made to this property.
    """
    PATH = "/workspace/flask/db"
    OUT_PATH = "/workspace/entap/outfiles"
    ENABLED_PATH = "/workspace/flask/databases.json"


    def __init__(
        self
    ):
        self.__files = listdir(self.PATH)
        with open(self.ENABLED_PATH,"r") as ifile:
            self.__enabled = loads(ifile.read())


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
                 index for all databases of this model, starting at 0 and ending
                 at the last database in the list.
        """
        return range(len(self.__files)).__iter__()


    def disable(
        self
        ,name
    ):
        """
        Disables the given database for use when running EnTAP. The model must
        be saved for this to be persistent.

        Parameters
        ----------
        name : string
               The file name of the database that is disabled.
        """
        if name in self.__enabled:
            self.__enabled.remove(name)


    def enable(
        self
        ,name
    ):
        """
        Enables the given database for use when running EnTAP. The model must be
        saved for this to be persistent.

        Parameters
        ----------
        name : string
               The file name of the database that is enabled.
        """
        if name not in self.__enabled:
            self.__enabled.append(name)


    def enabledNames(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : list
                 File names of all databases enabled in this model.
        """
        return self.__enabled


    def enabledPaths(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : list
                 Absolute paths of all databases enabled in this model's list.
        """
        return [self.__indexPath_(f) for f in self.__enabled]


    def extName(
        self
        ,index
    ):
        """
        Getter method.

        Parameters
        ----------
        index : int
                Index of the database whose extension name is returned.

        Returns
        -------
        result : string
                 Name of the extension of this model's database at the given
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
        Initializes this model's JSON configuration file that is used to store
        enabled properties of this model's databases. If the file does not exist
        it is created with an empty list of enabled databases.
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
        Getter method.

        Parameters
        ----------
        index : int
                Index of the database whose enabled state is returned.

        Returns
        -------
        result : bool
                 True if the given database is enabled or false otherwise.
        """
        return self.__files[index] in self.__enabled


    def isIndexed(
        self
        ,index
    ):
        """
        Getter method.

        Parameters
        ----------
        index : int
                Index of the database whose indexed state is returned.

        Returns
        -------
        result : bool
                 True if the given database is indexed or false otherwise.
        """
        return pathExists(self.__indexPath_(self.__files[index]))


    def name(
        self
        ,index
    ):
        """
        Getter method.

        Parameters
        ----------
        index : int
                Index of the database whose filename is returned.

        Returns
        -------
        result : string
                 Filename of this model's database at the given index.
        """
        return self.__files[index]


    def remove(
        self
        ,names
    ):
        """
        Removes the given list of databases from this model's list, including
        any diamond files of the given databases that are indexed. Saving this
        model is not required because the relevant database files themselves are
        removed.

        Parameters
        ----------
        names : list
                Filenames of databases that are removed from this model.
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
        Saves this model's enabled databases properties to its JSON file.
        """
        with open(self.ENABLED_PATH,"w") as ofile:
            ofile.write(dumps(self.__enabled))


    @classmethod
    def __indexPath_(
        cls
        ,name
    ):
        """
        Getter method.

        Parameters
        ----------
        name : string
               Filename of the database whose diamond index path is returned.

        Returns
        -------
        result : string
                 Full path of where the given database's diamond file would be
                 located. The diamond file may or may not exist.
        """
        return pathJoin(cls.OUT_PATH,"bin",name[:name.rfind(".")]+".dmnd")
