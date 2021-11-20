"""
Contains the RemoteUploadTask class.
"""
from ..Abstract.AbstractTask import *
from ..Model.DatabasesModel import *
from flask import render_template
from ftplib import FTP
import ftplib
from requests import get as reqGet
from os.path import exists as pathExists
from os.path import join as pathJoin




class RemoteUploadTask(AbstractTask):
    """
    This is the remote upload task class. It provides a task for this flask
    application to download databases from a list of provided URLs. It uses the
    path attribute from the databases model as the location where the databases
    are downloaded. The supported protocols are HTTP and FTP.
    """


    def __init__(
        self
        ,urls
    ):
        """
        Initializes this new remote upload task with the given list of URLs.

        Parameters
        ----------
        urls : list
               URLs of remote database locations that will be downloaded.
        """
        super().__init__()
        self.__urls = urls
        self.__url = ""
        self.__fileName = ""
        self.__path = ""
        self.__fi = 0
        self.__ft = len(urls)
        self.__finished = []
        self.__failed = []
        self._setRenderVars_(stage="init")


    def render(
        self
        ,**kwargs
    ):
        return render_template("task/remoteUpload.html",**kwargs)


    def run(
        self
    ):
        self.__finished = []
        self.__failed = []
        self.__fi = 1
        for url in self.__urls:
            self.__url = url
            self.__fileName = url[url.rfind("/")+1:]
            self._setRenderVars_(
                stage="prep"
                ,fileName=self.__fileName
                ,current=self.__fi
                ,total=self.__ft
            )
            try:
                self.__path = pathJoin(DatabasesModel.PATH,self.__fileName.replace(".*",""))
                if pathExists(self.__path):
                    self.__failed.append(
                        (
                            self.__fileName
                            ,"File already exists. Please remove it first if you wish to overwrite"
                                " it."
                        )
                    )
                    continue
                if url.lower().startswith("ftp://"):
                    self.__ftp_()
                else:
                    self.__query_()
                self.__fi += 1
            except Exception as e:
                self.__failed.append((self.__fileName,str(e)))
        self._setRenderVars_(stage="fin",finished=self.__finished,failed=self.__failed)
        return not self.__failed


    def title(
        self
    ):
        return "Remote Upload"


    def __ftp_(
        self
    ):
        """
        Downloads this task's currently active URL. The currently active URL
        must be a FTP protocol URL.
        """
        u = self.__url[6:]
        host = u[:u.find("/")]
        d = u[u.find("/"):u.rfind("/")+1]
        ftp = FTP(host)
        ftp.login()
        ftp.cwd(d)
        total = None
        if "*" not in self.__fileName:
            total = ftp.size(self.__fileName)
        if total is None:
            total = 0
        progress = 0
        with open(self.__path,"wb") as ofile:
            def write(chunk):
                nonlocal progress
                ofile.write(chunk)
                progress += len(chunk)
                p = (
                    int(min(100,progress*100/total)) if total
                    else self.__reportBySize_(progress)
                )
                self._setRenderVars_(
                    stage="down"
                    ,fileName=self.__fileName
                    ,current=self.__fi
                    ,total=self.__ft
                    ,progress=p
                    ,hasPercent=bool(total)
                )
            if "*" in self.__fileName:
                i = 1
                missed = 0
                while True:
                    try:
                        ftp.retrbinary("RETR "+self.__fileName.replace("*",str(i)),write)
                    except ftplib.error_perm as reason:
                        print(reason,missed)
                        if str(reason)[:3] != "550":
                            raise
                        else:
                            missed += 1
                            if missed > 10:
                                break
                    else:
                        missed = 0
                    finally:
                        i += 1
            else:
                ftp.retrbinary("RETR "+self.__fileName,write)
        self.__finished.append(self.__fileName)


    def __query_(
        self
    ):
        """
        Downloads this task's currently active URL. The currently active URL
        must be a HTTP protocol URL.
        """
        progress = 0
        rq = reqGet(self.__url,stream=True)
        total = int(rq.headers.get("Content-length",0))
        with open(self.__path,"wb") as ofile:
            for chunk in rq.iter_content(chunk_size=1024):
                if chunk:
                    ofile.write(chunk)
                    progress += 1024
                    p = (
                        int(min(100,progress*100/total)) if total
                        else self.__reportBySize_(progress)
                    )
                    self._setRenderVars_(
                        stage="down"
                        ,fileName=self.__fileName
                        ,current=self.__fi
                        ,total=self.__ft
                        ,progress=p
                        ,hasPercent=bool(total)
                    )
        self.__finished.append(self.__fileName)


    @staticmethod
    def __reportBySize_(
        size
    ):
        """
        Getter method.

        Parameters
        ----------
        size : int
               The raw size in bytes.

        Returns
        -------
        result : string
                 The size of the given number in proper units to keep the number
                 below 1024. These units are either bytes, kilobytes, megabytes,
                 gigabytes, or terabytes.
        """
        scale = ("B","KB","MB","GB","TB")
        i = 0
        while size > 1024 and i < len(scale):
            size /= 1024
            i += 1
        return f"{size:.2f}"+scale[i]
