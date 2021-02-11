"""
Contains the RemoteUploadTask class.
"""
from ..Abstract.AbstractTask import *
from ..Model.DatabasesModel import *
from flask import render_template
from ftplib import FTP
from requests import get as reqGet
from os.path import exists as pathExists
from os.path import join as pathJoin




class RemoteUploadTask(AbstractTask):
    """
    Detailed description.
    """


    def __init__(
        self
        ,urls
    ):
        """
        Detailed description.

        Parameters
        ----------
        urls : 
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
        """
        Detailed description.

        Parameters
        ----------
        **kwargs : 
        """
        return render_template("task/remoteUpload.html",**kwargs)


    def run(
        self
    ):
        """
        Detailed description.
        """
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
                self.__path = pathJoin(DatabasesModel.PATH,self.__fileName)
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
        """
        Detailed description.
        """
        return "Remote Upload"


    def __ftp_(
        self
    ):
        """
        Detailed description.
        """
        u = self.__url[6:]
        host = u[:u.find("/")]
        d = u[u.find("/"):u.rfind("/")+1]
        ftp = FTP(host)
        ftp.login()
        ftp.cwd(d)
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
            ftp.retrbinary("RETR "+self.__fileName,write)
        self.__finished.append(self.__fileName)


    def __query_(
        self
    ):
        """
        Detailed description.
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
        Detailed description.

        Parameters
        ----------
        size : 
        """
        scale = ("B","KB","MB","GB","TB")
        i = 0
        while size > 1024 and i < len(scale):
            size /= 1024
            i += 1
        return f"{size:.2f}"+scale[i]
