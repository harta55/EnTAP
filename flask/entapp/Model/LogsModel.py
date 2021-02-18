"""
Contains the LogsModel class.
"""
from .LogItem import *
from datetime import datetime
from os import listdir
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
                self.__logs.append(
                    LogItem(
                        self.PATH
                        ,datetime(*[int(n) for n in self.dateRE.findall(name)])
                    )
                )
        self.__logs.sort()
        self.__logs.reverse()


    def __iter__(
        self
    ):
        """
        Detailed description.
        """
        return range(len(self.__logs)).__iter__()


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
