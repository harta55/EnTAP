"""
Contains the DatabasesView class.
"""
from ..Controller import taskController
from ..Form.DatabaseUploadForm import *
from ..Model.DatabasesModel import *
from ..Task.RemoteUploadTask import *
import flask
from flask import flash
from flask import make_response
from flask import redirect
from flask import render_template
from flask import url_for
from flask_classful import FlaskView
from flask_classful import route
from os.path import join as pathJoin
from os.path import exists as pathExists
from os.path import getsize
from werkzeug.utils import secure_filename




class DatabasesView(FlaskView):
    """
    Detailed description.
    """
    route_base = "/databases/"


    def index(
        self
    ):
        """
        Detailed description.
        """
        databases = DatabasesModel()
        return render_template("databases/index.html",databases=databases)


    @route("/upload/chunk/",methods=["POST"])
    def uploadChunk(
        self
    ):
        """
        Detailed description.
        """
        ifile = flask.request.files["file"]
        savePath = pathJoin(DatabasesModel.PATH,secure_filename(ifile.filename))
        chunkIndex = int(flask.request.form["dzchunkindex"])
        if pathExists(savePath) and chunkIndex == 0:
            return make_response(("File already exists.",400))
        try:
            with open(savePath,"ab") as ofile:
                ofile.seek(int(flask.request.form["dzchunkbyteoffset"]))
                ofile.write(ifile.stream.read())
        except OSError:
            return make_response(("System error occured.",500))
        chunkTotal = int(flask.request.form["dztotalchunkcount"])
        if chunkIndex + 1 == chunkTotal:
            if getsize(savePath) != int(flask.request.form["dztotalfilesize"]):
                return make_response(("System error occured.",500))
        return make_response(("Chunk upload successful.",200))


    @route("/upload/local/")
    def uploadLocal(
        self
    ):
        """
        Detailed description.
        """
        return render_template("databases/uploadLocal.html")


    @route("/upload/remote/")
    def uploadRemote(
        self
    ):
        """
        Detailed description.
        """
        form = DatabaseUploadForm()
        return render_template("databases/uploadRemote.html",form=form)


    @route("/upload/remote/start/",methods=["POST"])
    def uploadRemoteStart(
        self
    ):
        """
        Detailed description.
        """
        form = DatabaseUploadForm()
        if form.validate():
            if form.submit.data:
                taskController.start(RemoteUploadTask(form.url.data))
                return redirect(url_for("RootView:index"))
        else:
            flash("Failed starting remote upload task.","danger")
        return redirect(url_for("DatabasesView:index"))
