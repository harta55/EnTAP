"""
Contains the GunzipTask class.
"""
from . import interfaces







class GunzipTask(interfaces.AbstractTask):
    """
    Detailed description.
    """


    def __init__(
        self
        ,workDir
        ,fileNames
        ):
        """
        Detailed description.

        Parameters
        ----------
        workDir : object
                  Detailed description.
        fileNames : object
                    Detailed description.
        """
        super().__init__()
        self.__result = enums.TaskResult.Running
        self.__workDir = workDir
        self.__names = fileNames
        self.__i = -1


    def errorOutput(
        self
        ):
        """
        Detailed description.
        """
        return self.__error


    def finalOutput(
        self
        ):
        """
        Detailed description.
        """
        return "Finished gunzipping all files."


    def output(
        self
        ):
        """
        Detailed description.
        """
        if self.__i >= 0:
            ret1 = "Gunzipping %i of %i\n" % (self.__i+1,len(self.__names))
            ret2 = self.__names[self.__i]
            return ret1+ret2


    def result(
        self
        ):
        """
        Detailed description.
        """
        return self.__result


    def run(
        self
        ):
        """
        Detailed description.
        """
        for name in self.__names:
            self.__i += 1
            if name.endswith(".gz"):
                path = os.path.join(self.__workDir,name)
                if os.path.isfile(path):
                    cmd = ["gunzip",path]
                    if subprocess.run(cmd).returncode != 0:
                        self.__error = "Failed gunzipping '" + name + "'."
                        self.__result = enums.TaskResult.Error
                        return
        self.__result = enums.TaskResult.Finished


    def title(
        self
        ):
        """
        Detailed description.
        """
        return "Gunzipping Files"
