"""
Contains the UninformsModel class.
"""
from json import dumps
from json import loads
from os import makedirs
from os.path import dirname
from os.path import exists as pathExists




class UninformsModel():
    """
    Detailed description.
    """
    PATH = "/workspace/flask/uninforms_config.json"


    def __init__(
        self
    ):
        """
        Detailed description.
        """
        with open(self.PATH,"r") as ifile:
            data = loads(ifile.read())
            self.__uninforms = set(data)


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
        return name in self.__uninforms


    def __iter__(
        self
    ):
        """
        Detailed description.
        """
        return self.__uninforms.__iter__()


    def add(
        self
        ,name
    ):
        """
        Detailed description.

        Parameters
        ----------
        name : 
        """
        if name in self.__uninforms:
            return False
        self.__uninforms.add(name)
        return True


    def configLines(
        self
    ):
        """
        Detailed description.
        """
        return ["uninformative="+u for u in self.__uninforms]


    @classmethod
    def initialize(
        cls
    ):
        """
        Detailed description.
        """
        if not pathExists(dirname(cls.PATH)):
            makedirs(dirname(cls.PATH))
        if not pathExists(cls.PATH):
            with open(cls.PATH,"w") as ofile:
                ofile.write(
                    dumps(
                        [
                            "conserved"
                            ,"predicted"
                            ,"unknown"
                            ,"hypothetical"
                            ,"putative"
                            ,"unidentified"
                            ,"uncultured"
                            ,"uninformative"
                        ]
                    )
                )


    def remove(
        self
        ,name
    ):
        """
        Detailed description.

        Parameters
        ----------
        name : 
        """
        self.__uninforms.remove(name)


    def save(
        self
    ):
        """
        Detailed description.
        """
        with open(self.PATH,"w") as ofile:
            ofile.write(dumps(list(self.__uninforms)))
