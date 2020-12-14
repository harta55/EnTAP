"""
Detailed description.
"""
from .. import core
import flask
from . import forms
from . import main
import os
from .. import tasks
from werkzeug import utils as wzutils




@main.route('/config/basic',methods=["GET","POST"])
def basicConfig():
    """
    Detailed description.
    """
    if core.taskManager.isRunning():
        return flask.redirect(flask.url_for("main.index"))
    form = forms.BasicConfigForm()
    if flask.request.method == "POST":
        if form.validate():
            form.save()
            flask.flash("Basic configuration successfully updated.","success")
            return flask.redirect(flask.url_for("main.basicConfig"))
    else:
        form.load()
    return flask.render_template('config/basic.html',form=form)




@main.route('/config/contams',methods=["GET","POST"])
def contamsConfig():
    """
    Detailed description.
    """
    if core.taskManager.isRunning():
        return flask.redirect(flask.url_for("main.index"))
    form = forms.ContamsConfigForm()
    if flask.request.method == "POST":
        if form.validate():
            form.processContams()
            form.save()
            return flask.redirect(flask.url_for("main.contamsConfig"))
    return flask.render_template('config/contams.html',form=form)




@main.route('/')
def index():
    """
    Detailed description.
    """
    core.taskManager.update()
    return flask.render_template('index.html',taskManager=core.taskManager)




@main.route('/run')
def run():
    """
    Detailed description.
    """
    if core.taskManager.isRunning():
        return flask.redirect(flask.url_for("main.index"))
    return flask.render_template('run.html')




@main.route('/config/uninforms',methods=["GET","POST"])
def uninformsConfig():
    """
    Detailed description.
    """
    if core.taskManager.isRunning():
        return flask.redirect(flask.url_for("main.index"))
    form = forms.UninformsConfigForm()
    if flask.request.method == "POST":
        if form.validate():
            form.processContams()
            form.save()
            return flask.redirect(flask.url_for("main.uninformsConfig"))
    return flask.render_template('config/uninforms.html',form=form)




def upload(
    workDir
    ):
    """
    Detailed description.

    Parameters
    ----------
    workDir : object
              Detailed description.
    """
    ifile = flask.request.files["file"]
    savePath = os.path.join(workDir,wzutils.secure_filename(ifile.filename))
    chunkIndex = int(flask.request.form["dzchunkindex"])
    if os.path.exists(savePath) and chunkIndex == 0:
        return flask.make_response(("File already exists.",400))
    try:
        with open(savePath,"ab") as ofile:
            ofile.seek(int(flask.request.form["dzchunkbyteoffset"]))
            ofile.write(ifile.stream.read())
    except OSError:
        return flask.make_response(("System error occured.",500))
    chunkTotal = int(flask.request.form["dztotalchunkcount"])
    if chunkIndex + 1 == chunkTotal:
        if os.path.getsize(savePath) != int(flask.request.form["dztotalfilesize"]):
            return flask.make_response(("System error occured.",500))
        else:
            #success
            pass
    return flask.make_response(("Chunk upload successful.",200))




@main.route('/upload/database',methods=["POST"])
def uploadDatabase():
    """
    Detailed description.
    """
    return upload("/workspace/flask/db")




@main.route('/upload/databases',methods=["GET","POST"])
def uploadDatabases():
    """
    Detailed description.
    """
    if core.taskManager.isRunning():
        return flask.redirect(flask.url_for("main.index"))
    remoteUploadForm = forms.RemoteDatabaseUploadForm()
    fileListForm = forms.FileListForm("/workspace/flask/db")
    if flask.request.method == "POST":
        if fileListForm.validate():
            fileListForm.removeSelected()
            return flask.redirect(flask.url_for("main.uploadDatabases"))
    return flask.render_template(
        "upload/database.html"
        ,fileListForm=fileListForm
        ,remoteUploadForm=remoteUploadForm
    )




@main.route('/upload/input',methods=["POST"])
def uploadInput():
    """
    Detailed description.
    """
    return upload("/workspace/entap/entap_infiles")




@main.route('/upload/inputs',methods=["GET","POST"])
def uploadInputs():
    """
    Detailed description.
    """
    if core.taskManager.isRunning():
        return flask.redirect(flask.url_for("main.index"))
    form = forms.FileListForm("/workspace/entap/entap_infiles")
    if flask.request.method == "POST":
        if form.validate():
            form.removeSelected()
            return flask.redirect(flask.url_for("main.uploadInputs"))
    return flask.render_template("upload/input.html",form=form)




@main.route('/upload/remote/database',methods=["POST"])
def uploadRemoteDatabase():
    """
    Detailed description.
    """
    if core.taskManager.isRunning():
        return flask.redirect(flask.url_for("main.index"))
    form = forms.RemoteDatabaseUploadForm()
    if form.validate():
        print(form.customURL.data)
        task = tasks.DownloadTask(form.customURL.data)
        core.taskManager.start(task)
    return flask.redirect(flask.url_for("main.uploadDatabases"))