"""
Contains the RemoteUploadTask class.
"""
from ..Abstract.AbstractTask import *
from ..Model.DatabasesModel import *
from flask import render_template
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
        self._setRenderVars_(running=True,fileName="",current=0,total=0,progress=0,totalSize=True)


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
        finished = []
        failed = []
        fi = 1
        ft = len(self.__urls)
        for url in self.__urls:
            try:
                fileName = url[url.rfind("/")+1:]
                progress = 0
                rq = reqGet(url,stream=True)
                total = int(rq.headers.get("Content-length",0))
                path = pathJoin(DatabasesModel.PATH,fileName)
                if pathExists(path):
                    failed.append(
                        (
                            fileName
                            ,"File already exists. Please remove it first if you wish to overwrite"
                             " it."
                        )
                    )
                    continue
                with open(path,"wb") as ofile:
                    for chunk in rq.iter_content(chunk_size=1024):
                        if chunk:
                            ofile.write(chunk)
                            progress += 1024
                            p = (
                                int(min(100,progress*100/total)) if total
                                else self.__reportBySize_(progress)
                            )
                            self._setRenderVars_(
                                running=True
                                ,fileName=fileName
                                ,current=fi
                                ,total=ft
                                ,progress=p
                                ,totalSize=bool(total)
                            )
                fi += 1
                finished.append(fileName)
            except Exception as e:
                failed.append((fileName,str(e)))
        self._setRenderVars_(running=False,finished=finished,failed=failed)
        return not failed


    def title(
        self
    ):
        """
        Detailed description.
        """
        return "Remote Upload"


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
