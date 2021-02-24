"""
Contains the LogItem class.
"""
from os.path import join as pathJoin




class LogItem():
    """
    This is the log item class. It provides a model for a single EnTAP log item.
    It is a read only model that provides the date, log, debug log, and tail
    preview of the log.
    """


    def __init__(
        self
        ,workDir
        ,dt
    ):
        """
        Initializes this new log item with the given working directory and date.

        Parameters
        ----------
        workDir : string
                  Path of the working directory where all EnTAP log files are
                  located.
        dt : datetime.datetime
             The date and time for this item's EnTAP log.
        """
        self.__workDir = workDir
        self.__dt = dt


    def __lt__(
        self
        ,other
    ):
        return self.__dt.__lt__(other.__dt)


    def debug(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : string
                 Contents of this item's EnTAP debug log file.
        """
        x = self.__dt
        path = pathJoin(
            self.__workDir
            ,"debug_%iY%iM%iD-%ih%im%is.txt"%(x.year,x.month,x.day,x.hour,x.minute,x.second)
        )
        with open(path,"r") as ifile:
            return ifile.read()


    def dt(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : datetime.datetime
                 The date and time of this item's EnTAP log.
        """
        return self.__dt


    def read(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : string
                 Contents of this item's EnTAP log file.
        """
        x = self.__dt
        path = pathJoin(
            self.__workDir
            ,"log_file_%iY%iM%iD-%ih%im%is.txt"%(x.year,x.month,x.day,x.hour,x.minute,x.second)
        )
        with open(path,"r") as ifile:
            return ifile.read()


    def tail(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : string
                 The last ten lines of this item's EnTAP log file.
        """
        x = self.__dt
        path = pathJoin(
            self.__workDir
            ,"log_file_%iY%iM%iD-%ih%im%is.txt"%(x.year,x.month,x.day,x.hour,x.minute,x.second)
        )
        with open(path,"r") as ifile:
            lines = ifile.read().split("\n")
            return "\n".join(lines[-10:])
