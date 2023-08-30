"""
Contains the ContamsModel class.
"""
from json import dumps
from json import loads
from os import makedirs
from os.path import dirname
from os.path import exists as pathExists




class ContamsModel():
    """
    This is the contaminants model class. It provides a model for the
    contaminants list configuration of EnTAP. It interfaces with its respective
    form class to populate the form or update values from the form. This class
    can generate its portion of EnTAP INI configuration output.
    """
    PATH = "/workspace/flask/contams_config.json"


    def __init__(
        self
    ):
        with open(self.PATH,"r") as ifile:
            data = loads(ifile.read())
            self.__contams = set(data)


    def __contains__(
        self
        ,name
    ):
        return name in self.__contams


    def __iter__(
        self
    ):
        return self.__contams.__iter__()


    def add(
        self
        ,name
    ):
        """
        Adds a new contaminant to this model's list with the given name.

        Parameters
        ----------
        name : string
               Name of the new contaminant.
        """
        if name in self.__contams:
            return False
        self.__contams.add(name)
        return True


    def configLines(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : list
                 Contaminants configuration for the EnTAP INI file generated
                 from the current settings of this model. Each string represents
                 a line of output.
        """
        return ["contam=" + ",".join([c.replace(" ","_") for c in self.__contams])]


    @classmethod
    def initialize(
        cls
    ):
        """
        Initializes this model's JSON configuration file that is used to store
        the values of this configuration model. If the file does not exist it is
        created with default values.
        """
        if not pathExists(dirname(cls.PATH)):
            makedirs(dirname(cls.PATH))
        if not pathExists(cls.PATH):
            with open(cls.PATH,"w") as ofile:
                ofile.write(dumps([]))


    def remove(
        self
        ,name
    ):
        """
        Removes the given contaminant name from this model's list. The given
        contaminant name must exist in this model's list.

        Parameters
        ----------
        name : string
               Name of the contaminant that is removed.
        """
        self.__contams.remove(name)


    def save(
        self
    ):
        """
        Saves this model's current values to its JSON file.
        """
        with open(self.PATH,"w") as ofile:
            ofile.write(dumps(list(self.__contams)))
