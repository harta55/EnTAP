"""
Detailed description.
"""
import flask
import flask_bootstrap




def create_app():
    """
    Detailed description.
    """
    app = flask.Flask(__name__)
    flask_bootstrap.Bootstrap(app)
    from .main import main as mainBlueprint
    app.register_blueprint(mainBlueprint)
    return app
