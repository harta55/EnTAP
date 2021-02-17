"""
Contains the LogItem class.
"""




class LogItem():
    """
    Detailed description.
    """


    def __init__(
        self
        ,workDir
        ,dt
    ):
        """
        Detailed description.

        Parameters
        ----------
        workDir : 
        dt : 
        """
        self.__workDir = workDir
        self.__dt = dt


    def __lt__(
        self
        ,other
    ):
        """
        Detailed description.

        Parameters
        ----------
        other : 
        """
        return self.__dt.__lt__(other.__dt)


    def debug(
        self
    ):
        """
        Detailed description.
        """
        path = pathJoin(
            self.__workDir
            ,"debug_%iY%iM%iD-%id%ih%im.txt"%(x.year,x.month,x.day,x.hour,x.minute,x.second)
        )
        with open(path,"r") as ifile:
            return ifile.read()


    def dt(
        self
    ):
        """
        Detailed description.
        """
        return self.__dt


    def read(
        self
    ):
        """
        Detailed description.
        """
        path = pathJoin(
            self.__workDir
            ,"log_file_%iY%iM%iD-%id%ih%im.txt"%(x.year,x.month,x.day,x.hour,x.minute,x.second)
        )
        with open(path,"r") as ifile:
            return ifile.read()


    def tail(
        self
    ):
        """
        Detailed description.
        """
        path = pathJoin(
            self.__workDir
            ,"log_file_%iY%iM%iD-%id%ih%im.txt"%(x.year,x.month,x.day,x.hour,x.minute,x.second)
        )
        with open(path,"r") as ifile:
            lines = ifile.read().split("\n")
            return "\n".join(lines[-10:])
