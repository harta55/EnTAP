"""
Detailed description.
"""
import flask
from . import forms
from . import main




@main.route('/config/basic',methods=["GET","POST"])
def basicConfig():
    """
    Detailed description.
    """
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
    return flask.render_template('index.html')




@main.route('/run')
def run():
    """
    Detailed description.
    """
    return flask.render_template('run.html')




@main.route('/config/uninforms',methods=["GET","POST"])
def uninformsConfig():
    """
    Detailed description.
    """
    form = forms.UninformsConfigForm()
    if flask.request.method == "POST":
        if form.validate():
            form.processContams()
            form.save()
            return flask.redirect(flask.url_for("main.uninformsConfig"))
    return flask.render_template('config/uninforms.html',form=form)




@main.route('/upload/database')
def uploadDatabases():
    """
    Detailed description.
    """
    return flask.render_template('upload/database.html')




@main.route('/upload/input')
def uploadInputs():
    """
    Detailed description.
    """
    return flask.render_template('upload/input.html')
