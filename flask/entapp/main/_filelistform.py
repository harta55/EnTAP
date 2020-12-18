"""
Contains the FileListForm class.
"""
import flask
import flask_wtf
import os
import wtforms as wtf








class FileListForm(flask_wtf.FlaskForm):
    """
    Detailed description.
    """
    fileList = wtf.SelectMultipleField(
        "File List"
        ,choices=[]
        ,description="Select any number of listed files."
    )
    gzipSubmit = wtf.SubmitField("Gunzip Selected Files")
    rmSubmit = wtf.SubmitField("Remove Selected Files")


    def __init__(
        self
        ,workDir
        ):
        """
        Detailed description.

        Parameters
        ----------
        workDir : object
                  Detailed description.
        """
        super().__init__()
        self.__workDir = workDir
        self.fileList.choices = [(f,f) for f in os.listdir(workDir)]


    def processAll(
        self
        ):
        """
        Detailed description.
        """
        if self.rmSubmit:
            self.__removeSelected
        elif self.gzipSubmit:
            self.__gunzipSelected


    def __gunzipSelected(
        self
        ):
        """
        Detailed description.
        """
        pass


    def __removeSelected(
        self
        ):
        """
        Detailed description.
        """
        if self.fileList.data:
            for path in self.fileList.data:
                path = os.path.join(self.__workDir,path)
                if os.path.isfile(path):
                    os.remove(path)
            flask.flash("Selected files successfully removed.","success")
