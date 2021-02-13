"""
Contains the DatabasesView class.
"""
from ..Controller import taskController
from ..Form.DatabaseUploadForm import *
from ..Model.DatabasesModel import *
from ..Task.IndexTask import *
from ..Task.RemoteUploadTask import *
from flask import flash
from flask import make_response
from flask import redirect
from flask import render_template
from flask import request
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


    @route("/operation/",methods=["POST"])
    def operation(
        self
    ):
        """
        Detailed description.
        """
        action = request.form.get("action")
        names = request.form.getlist("names")
        if not names:
            return redirect(url_for("DatabasesView:index"))
        if action == "index":
            taskController.start(IndexTask(names))
            flash("Started index task.","success")
            return redirect(url_for("RootView:status"))
        elif action == "enable":
            databases = DatabasesModel()
            for name in names:
                databases.enable(name)
            databases.save()
            flash("Successfully enabled database(s).","success")
            return redirect(url_for("DatabasesView:index"))
        elif action == "disable":
            databases = DatabasesModel()
            for name in names:
                databases.disable(name)
            databases.save()
            flash("Successfully disabled database(s).","success")
            return redirect(url_for("DatabasesView:index"))
        elif action == "remove":
            databases = DatabasesModel()
            if databases.remove(names):
                flash("Successfully removed database(s).","success")
            else:
                flash("Failed removing one or more databases.","danger")
            return redirect(url_for("DatabasesView:index"))
        else:
            flash("Unknown operation given.","danger")
            return redirect(url_for("DatabasesView:index"))


    @route("/upload/chunk/",methods=["POST"])
    def uploadChunk(
        self
    ):
        """
        Detailed description.
        """
        ifile = request.files["file"]
        savePath = pathJoin(DatabasesModel.PATH,secure_filename(ifile.filename))
        chunkIndex = int(request.form["dzchunkindex"])
        if pathExists(savePath) and chunkIndex == 0:
            return make_response(("File already exists.",400))
        try:
            with open(savePath,"ab") as ofile:
                ofile.seek(int(request.form["dzchunkbyteoffset"]))
                ofile.write(ifile.stream.read())
        except OSError:
            return make_response(("System error occured.",500))
        chunkTotal = int(request.form["dztotalchunkcount"])
        if chunkIndex + 1 == chunkTotal:
            if getsize(savePath) != int(request.form["dztotalfilesize"]):
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
                urls = [u.strip() for u in form.urls.data.split("\n") if u]
                urls += form.defUrls.data
                taskController.start(RemoteUploadTask(urls))
                flash("Started remote upload task.","success")
                return redirect(url_for("RootView:status"))
        else:
            flash("Failed starting remote upload task.","danger")
        return redirect(url_for("DatabasesView:index"))
