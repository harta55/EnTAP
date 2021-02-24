"""
Contains the LogsModel class.
"""
from .LogItem import *
from datetime import datetime
from os import listdir
from re import compile as reCompile




class LogsModel():
    """
    This is the logs model class. It provides a read only model for all EnTAP
    logs as a list of log item class instances.
    """
    PATH = "/workspace/entap/outfiles"
    dateRE = reCompile("\d+")


    def __init__(
        self
    ):
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
        return range(len(self.__logs)).__iter__()


    def __getitem__(
        self
        ,index
    ):
        return self.__logs[index]


    def __len__(
        self
    ):
        return len(self.__logs)
