"""
Contains the AbstractTask class.
"""
import abc
import threading








class AbstractTask(abc.ABC, threading.Thread):
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


    @abc.abstractmethod
    def errorOutput(
        self
        ):
        """
        Detailed description.
        """
        pass


    @abc.abstractmethod
    def finalOutput(
        self
        ):
        """
        Detailed description.
        """
        pass


    @abc.abstractmethod
    def output(
        self
        ):
        """
        Detailed description.
        """
        pass


    @abc.abstractmethod
    def result(
        self
        ):
        """
        Detailed description.
        """
        pass


    @abc.abstractmethod
    def title(
        self
        ):
        """
        Detailed description.
        """
        pass
