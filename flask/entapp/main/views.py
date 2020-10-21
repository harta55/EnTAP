"""
Detailed description.
"""
import flask
from . import forms
from . import main




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




@main.route('/config',methods=["GET","POST"])
def config():
    """
    Detailed description.
    """
    form = forms.BasicConfigForm()
    if flask.request.method == "POST":
        if form.validate():
            form.save()
            flask.flash("success","Configuration successfully updated.")
            return flask.redirect(flask.url_for("main.config"))
    else:
        form.load()
    return flask.render_template('config.html',form=form)
