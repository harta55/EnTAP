"""
Detailed description.
"""
import flask
import json
import os

from ._taskmanager import TaskManager




def create_app():
    """
    Detailed description.
    """
    app = flask.Flask(__name__)
    app.config["SECRET_KEY"] = "hard to guess phrase DO NOT USE IN PRODUCTION"
    from .main import main as mainBlueprint
    app.register_blueprint(mainBlueprint)
    if not os.path.exists("/workspace/flask"):
        os.makedirs("/workspace/flask/db")
        os.makedirs("/workspace/entap/entap_infiles")
        with open("/workspace/flask/basic_config.json","w") as ofile:
            ofile.write(
                json.dumps(
                    {
                        "dataType": "1"
                        ,"fpkm": "0.50"
                        ,"singleEnd": True
                        ,"seqComplete": True
                        ,"minProteinLength": 100
                        ,"refineStarts": True
                        ,"outputFormat": [
                            "1"
                            ,"3"
                            ,"4"
                        ]
                        ,"taxonomy": ""
                        ,"queryCoverage": "50.00"
                        ,"targetCoverage": "50.00"
                        ,"eValue": "0.000001"
                    }
                        ,indent=4
                )
            )
        with open("/workspace/flask/contams_config.json","w") as ofile:
            ofile.write(json.dumps([],indent=4))
        with open("/workspace/flask/uninforms_config.json","w") as ofile:
            ofile.write(json.dumps([],indent=4))
    return app
