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
        self.__output = []
        self.__finalOutput = ""
        self.__lock = Lock()


    def output(
        self
    ):
        """
        Detailed description.
        """
        self.__lock.acquire()
        ret = self.__finalOutput
        self.__lock.release()
        return ret


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


    def _clear_(
        self
    ):
        """
        Detailed description.
        """
        self.__output = []


    def _print_(
        self
        ,text
    ):
        """
        Detailed description.

        Parameters
        ----------
        text : 
        """
        self.__output.append(text)


    def _flush_(
        self
    ):
        """
        Detailed description.
        """
        self.__lock.acquire()
        self.__finalOutput = "".join(self.__output)
        self.__lock.release()
