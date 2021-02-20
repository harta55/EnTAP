"""
Contains the TaskController class.
"""
from ..Abstract.AbstractTask import *
from threading import Thread
from traceback import print_exc




class TaskController():
    """
    This is the singleton task controller. It manages tasks for this flask
    application, starting new tasks and providing information about the current
    state of any running task.
    """


    def __init__(
        self
    ):
        self.__thread = None
        self.__task = None
        self.__lastTitle = "No History"
        self.__lastOutput = "No task has been started."
        self.__hadError = False


    def hadError(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : bool
                 True if the last task ran failed or false otherwise. If no task
                 has been ran yet then false is returned.
        """
        return self.__hadError


    def isRunning(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : bool
                 True if this controller is running a task or false otherwise.
        """
        return self.__thread and self.__thread.is_alive()


    def output(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : string
                 The rendered output of this controller's currently running task
                 or the last task ran if no task is running. If no task has ever
                 been run a default output stating so is returned.
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
        Starts a new task for this controller in its separate task thread. This
        controller must not be running a task.

        Parameters
        ----------
        task : entapp.Abstract.AbstractTask
               New task that this controller will start running.
        """
        assert(not self.isRunning())
        assert(isinstance(task,AbstractTask))
        self.__hadError = False
        self.__thread = Thread(target=self.__run_)
        self.__task = task
        self.__thread.start()


    def title(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : string
                 The title of this controller's currently running task or the
                 last task ran if no task is running. If no task has ever been
                 run a default title stating so is returned.
        """
        if self.isRunning():
            return self.__task.title()
        else:
            return self.__lastTitle


    def update(
        self
    ):
        """
        Updates this task controller, checking if any currently running task is
        finished.
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
        Runs this controller's current task on its separate task thread,
        catching any thrown exception and marking the task had an error.
        """
        try:
            self.__hadError = not self.__task.run()
        except Exception:
            print_exc()
            self.__hadError = True
