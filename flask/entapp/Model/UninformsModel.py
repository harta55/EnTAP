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
    This is the uninformative model class. It provides a model for the
    uninformative list configuration of EnTAP. It interfaces with its respective
    form class to populate the form or update values from the form. This class
    can generate its portion of EnTAP INI configuration output.
    """
    PATH = "/workspace/flask/uninforms_config.json"


    def __init__(
        self
    ):
        with open(self.PATH,"r") as ifile:
            data = loads(ifile.read())
            self.__uninforms = set(data)


    def __contains__(
        self
        ,name
    ):
        return name in self.__uninforms


    def __iter__(
        self
    ):
        return self.__uninforms.__iter__()


    def add(
        self
        ,name
    ):
        """
        Adds a new uninformative to this model's list with the given name.

        Parameters
        ----------
        name : string
               Name of the new uninformative.
        """
        if name in self.__uninforms:
            return False
        self.__uninforms.add(name)
        return True


    def configLines(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : list
                 Uninformative configuration for the EnTAP INI file generated
                 from the current settings of this model. Each string represents
                 a line of output.
        """
        return ["uninformative="+u for u in self.__uninforms]


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
        Removes the given uninformative name from this model's list. The given
        uninformative name must exist in this model's list.

        Parameters
        ----------
        name : string
               Name of the uninformative that is removed.
        """
        self.__uninforms.remove(name)


    def save(
        self
    ):
        """
        Saves this model's current values to its JSON file.
        """
        with open(self.PATH,"w") as ofile:
            ofile.write(dumps(list(self.__uninforms)))
