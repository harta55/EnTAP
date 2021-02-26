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
    This is the databases view. It provides a flask view for the databases
    interface of this flask application. Because EnTAP databases are complex
    this class provides more than a few trivial routes. Services provided
    include local uploads, remote uploads, and various operations on contained
    databases.
    """
    route_base = "/databases/"


    def index(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : object
                 Index page of this view.
        """
        databases = DatabasesModel()
        return render_template("databases/index.html",databases=databases)


    @route("/operation/",methods=["POST"])
    def operation(
        self
    ):
        """
        Runs an operation on one or more databases this flask application
        contains in its local container. The action and list of databases to do
        it on is provided by the submitted form. The possible actions are index,
        enable, disable, and remove.

        Returns
        -------
        result : object
                 Flask redirect to the root view's status page if the action is
                 to index otherwise this view's index page.
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
        Processes a single chunk upload of a database file using the local
        upload interface.

        Returns
        -------
        result : object
                 Flask response of success or failure.
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
        Getter method.

        Returns
        -------
        result : object
                 This view's local upload page.
        """
        return render_template("databases/uploadLocal.html")


    @route("/upload/remote/")
    def uploadRemote(
        self
    ):
        """
        Getter method.

        Returns
        -------
        result : object
                 This view's remote upload page.
        """
        form = DatabaseUploadForm()
        return render_template("databases/uploadRemote.html",form=form)


    @route("/upload/remote/start/",methods=["POST"])
    def uploadRemoteStart(
        self
    ):
        """
        Starts a new remote upload task with the submitted form's list of URLs
        is it passes validation.

        Returns
        -------
        result : object
                 Flask redirect to the root view's status page if the remote
                 upload was successfully started otherwise this view's index
                 page.
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
