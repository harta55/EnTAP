"""
Contains the TaskManager class.
"""
from . import enums








class TaskManager():
    """
    Detailed description.
    """


    def __init__(
        self
        ):
        """
        Detailed description.
        """
        self.__task = None
        self.__result = enums.TaskResult.Finished


    def isRunning(
        self
        ):
        """
        Detailed description.
        """
        return self.state() == enums.TaskManagerState.Running


    def output(
        self
        ):
        """
        Detailed description.
        """
        if self.isRunning():
            return self.__task.output()
        else:
            return (None,None)


    def start(
        self
        ,task
        ):
        """
        Detailed description.

        Parameters
        ----------
        task : object
               Detailed description.
        """
        if not self.isRunning():
            self.__task = task
            self.__task.start()


    def state(
        self
        ):
        """
        Detailed description.
        """
        if self.__task is None:
            if self.__result == enums.TaskResult.Finished:
                return enums.TaskManagerState.Idle
            else:
                return enums.TaskManagerState.Error
        else:
            return enums.TaskManagerState.Running


    def update(
        self
        ):
        """
        Detailed description.
        """
        if self.isRunning():
            if not self.__task.is_alive():
                self.__result = self.__task.result()
                self.__task = None
