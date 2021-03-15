"""
Contains the EnTAP flask application package. This package creates the flask app
in its root namespace with the name 'application'.
"""
from .Application import *


application = Application(__name__)
application.config["SECRET_KEY"] = "hard to guess phrase DO NOT USE IN PRODUCTION"


from .View.ConfigView import *
from .View.ContamsView import *
from .View.DatabasesView import *
from .View.InputsView import *
from .View.LogsView import *
from .View.RootView import *
from .View.RunView import *
from .View.UninformsView import *


ConfigView.register(application)
ContamsView.register(application)
DatabasesView.register(application)
InputsView.register(application)
LogsView.register(application)
RootView.register(application)
RunView.register(application)
UninformsView.register(application)


from . import Controller


@application.before_request
def update():
    ret = Controller.update()
    if ret:
        return ret
