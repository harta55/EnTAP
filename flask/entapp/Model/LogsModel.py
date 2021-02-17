"""
Contains the LogsModel class.
"""
from .LogItem import *
from os import listdir
from os.path import join as pathJoin
from re import compile as reCompile




class LogsModel():
    """
    Detailed description.
    """
    PATH = "/workspace/entap/outfiles"
    dateRE = reCompile("\d+")


    def __init__(
        self
    ):
        """
        Detailed description.
        """
        self.__logs = []
        for name in listdir(self.PATH):
            if name.startswith("log_file_") and name.endswith(".txt"):
                self.__logs.append(LogItem(pathJoin(self.PATH,datetime(*dateRE.findall(name)))))
        self.__logs.sort()


    def __iter__(
        self
    ):
        """
        Detailed description.
        """
        return self.__logs.__iter__()


    def __getitem__(
        self
        ,index
    ):
        """
        Detailed description.

        Parameters
        ----------
        index : 
        """
        return self.__logs[index]


    def __len__(
        self
    ):
        """
        Detailed description.
        """
        return len(self.__logs)
