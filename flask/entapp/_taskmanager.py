"""
Contains the TaskManager class.
"""
from . import enums
from . import interfaces








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
        self.__lastTask = None
        self.__result = enums.TaskResult.Finished


    def hasError(
        self
        ):
        """
        Detailed description.
        """
        return self.state() == enums.TaskManagerState.Error


    def isFinished(
        self
        ):
        """
        Detailed description.
        """
        return self.state() == enums.TaskManagerState.Idle and self.__lastTask is not None


    def isIdle(
        self
        ):
        """
        Detailed description.
        """
        return self.state() == enums.TaskManagerState.Idle and self.__lastTask is None


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
            return (self.__task.title(),self.__task.output())
        elif self.isFinished():
            return (self.__lastTask.title(),self.__lastTask.finalOutput())
        elif self.hasError():
            return (self.__lastTask.title(),self.__lastTask.errorOutput())
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
        assert(isinstance(task,interfaces.AbstractTask))
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
                self.__lastTask = self.__task
                self.__task = None
