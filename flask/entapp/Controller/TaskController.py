"""
Contains the TaskController class.
"""
from ..Abstract.AbstractTask import *
from threading import Thread




class TaskController():
    """
    Detailed description.
    """


    def __init__(
        self
    ):
        """
        Detailed description.
        """
        self.__thread = None
        self.__task = None
        self.__lastTitle = "Idle"
        self.__lastOutput = "No task has been started."
        self.__hadError = False


    def hadError(
        self
    ):
        """
        Detailed description.
        """
        return self.__hadError


    def isRunning(
        self
    ):
        """
        Detailed description.
        """
        return self.__thread and self.__thread.is_alive()


    def output(
        self
    ):
        """
        Detailed description.
        """
        if self.isRunning():
            return self.__task.output()
        else:
            return self.__lastOutput


    def start(
        self
        ,task
    ):
        """
        Detailed description.

        Parameters
        ----------
        task : 
        """
        assert(not self.isRunning())
        assert(isinstance(task,AbstractTask))
        self.__thread = Thread(target=self.__run_)
        self.__task = task
        self.__thread.start()


    def title(
        self
    ):
        """
        Detailed description.
        """
        if self.isRunning():
            return self.__task.title()
        else:
            return self.__lastTitle


    def update(
        self
    ):
        """
        Detailed description.
        """
        if self.__thread is not None:
            if not self.__thread.is_alive():
                self.__thread = None
                self.__lastTitle = self.__task.title()
                self.__lastOutput = self.__task.output()
                self.__task = None


    def __run_(
        self
    ):
        """
        Detailed description.
        """
        self.__hadError = not self.__task.run()
