"""
Detailed description.
"""
import flask




def create_app():
    """
    Detailed description.
    """
    app = flask.Flask(__name__)
    app.config["SECRET_KEY"] = "hard to guess phrase DO NOT USE IN PRODUCTION"
    from .main import main as mainBlueprint
    app.register_blueprint(mainBlueprint)
    return app
