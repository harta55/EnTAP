"""
Contains the InputsView class.
"""
from ..Controller import taskController
from ..Model.InputsModel import *
from ..Task.GunzipTask import *
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




class InputsView(FlaskView):
    """
    Detailed description.
    """
    route_base = "/inputs/"


    def index(
        self
    ):
        """
        Detailed description.
        """
        inputs = InputsModel()
        return render_template("inputs/index.html",inputs=inputs)


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
            return redirect(url_for("InputsView:index"))
        if action == "remove":
            inputs = InputsModel()
            if inputs.remove(names):
                flash("Successfully removed inputs.","success")
            else:
                flash("Failed removing one or more inputs.","danger")
            return redirect(url_for("InputsView:index"))
        elif action == "gunzip":
            taskController.start(GunzipTask(names))
            flash("Started gunzip task.","success")
            return redirect(url_for("RootView:status"))
        else:
            flash("Unknown operation given.","danger")
            return redirect(url_for("InputsView:index"))


    @route("/upload/chunk/",methods=["POST"])
    def uploadChunk(
        self
    ):
        """
        Detailed description.
        """
        ifile = request.files["file"]
        savePath = pathJoin(InputsModel.PATH,secure_filename(ifile.filename))
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
        return render_template("inputs/uploadLocal.html")
