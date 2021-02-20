"""
Contains the AbstractTask class.
"""
from abc import ABC
from abc import abstractmethod
from threading import Lock




class AbstractTask(ABC):
    """
    This is the abstract task class. It is a process interface designed to run a
    specific task in the background of this flask application on a separate
    thread to prevent locking up the main thread providing server responses.
    Because flask's render template only works on the main server thread, one
    interface is provided for rendering templates which is called from the main
    server thread. Variables can be passed from the task thread to the server
    thread with another provided method.
    """


    def __init__(
        self
    ):
        super().__init__()
        self.__renderVars = {}
        self.__lock = Lock()


    def output(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : string
                 Rendered HTML that is displayed to the user showing them the
                 status of this task.
        """
        self.__lock.acquire()
        ret = self.render(**self.__renderVars)
        self.__lock.release()
        return ret


    @abstractmethod
    def render(
        self
        ,**kwargs
    ):
        """
        This interface is a getter method. This is called in the main server
        thread.

        Parameters
        ----------
        **kwargs : dictionary
                   Keyword arguments last set by this task using the set render
                   variables method.

        Returns
        -------
        result : string
                 Rendered HTML that is displayed to the user showing them the
                 status of this task.
        """
        pass


    @abstractmethod
    def run(
        self
    ):
        """
        This interface runs this task. This is called in a separate task thread.
        """
        pass


    @abstractmethod
    def title(
        self
    ):
        """
        This interface is a getter method.

        Returns
        -------
        result : string
                 The title of this task.
        """
        pass


    def _setRenderVars_(
        self
        ,**kwargs
    ):
        """
        Sets the keyword render variables that are given to the render interface
        whenever it is called from the main server thread.

        Parameters
        ----------
        **kwargs : dictionary
                   New keyword arguments given to the render interface whenever
                   it is called.
        """
        self.__lock.acquire()
        self.__renderVars = kwargs
        self.__lock.release()
