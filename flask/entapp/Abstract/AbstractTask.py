"""
Contains the AbstractTask class.
"""
from abc import ABC
from abc import abstractmethod
from threading import Lock




class AbstractTask(ABC):
    """
    Detailed description.
    """


    def __init__(
        self
    ):
        """
        Detailed description.
        """
        super().__init__()
        self.__renderVars = {}
        self.__lock = Lock()


    def output(
        self
    ):
        """
        Detailed description.
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
        Detailed description.

        Parameters
        ----------
        **kwargs : 
        """
        pass


    @abstractmethod
    def run(
        self
    ):
        """
        Detailed description.
        """
        pass


    @abstractmethod
    def title(
        self
    ):
        """
        Detailed description.
        """
        pass


    def _setRenderVars_(
        self
        ,**kwargs
    ):
        """
        Detailed description.

        Parameters
        ----------
        **kwargs : 
        """
        self.__lock.acquire()
        self.__renderVars = kwargs
        self.__lock.release()
