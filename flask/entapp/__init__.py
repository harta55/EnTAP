"""
Detailed description.
"""
from .Application import *


application = Application(__name__)
application.config["SECRET_KEY"] = "hard to guess phrase DO NOT USE IN PRODUCTION"


from .View.ConfigView import *
from .View.ContamsView import *
from .View.DatabasesView import *
from .View.RootView import *
from .View.UninformsView import *


ConfigView.register(application)
ContamsView.register(application)
DatabasesView.register(application)
RootView.register(application)
UninformsView.register(application)


from . import Controller


@application.before_request
def update():
    ret = Controller.update()
    if ret:
        return ret
