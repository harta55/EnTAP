"""
Detailed description.
"""
import flask
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




@main.route('/config')
def config():
    """
    Detailed description.
    """
    return flask.render_template('config.html')
