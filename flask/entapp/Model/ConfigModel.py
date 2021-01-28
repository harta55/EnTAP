"""
Contains the ConfigModel class.
"""
from decimal import Decimal
from json import dumps
from json import loads
from os import makedirs
from os.path import dirname
from os.path import exists as pathExists




class ConfigModel():
    """
    Detailed description.
    """
    PATH = "/workspace/flask/config.json"


    def __init__(
        self
    ):
        """
        Detailed description.
        """
        with open(self.PATH,"r") as ifile:
            self.__data = loads(ifile.read())


    @classmethod
    def initialize(
        cls
    ):
        """
        Detailed description.
        """
        if not pathExists(dirname(cls.PATH)):
            makedirs(dirname(cls.PATH))
        if not pathExists(cls.PATH):
            with open(cls.PATH,"w") as ofile:
                ofile.write(
                    dumps(
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


    def populate(
        self
        ,form
    ):
        """
        Detailed description.

        Parameters
        ----------
        form : 
        """
        form.dataType.data = self.__data["dataType"]
        form.fpkm.data = Decimal(self.__data["fpkm"])
        form.singleEnd.data = self.__data["singleEnd"]
        form.seqComplete.data = self.__data["seqComplete"]
        form.minProteinLength.data = self.__data["minProteinLength"]
        form.refineStarts.data = self.__data["refineStarts"]
        form.outputFormat.data = self.__data["outputFormat"]
        form.taxonomy.data = self.__data["taxonomy"]
        form.queryCoverage.data = Decimal(self.__data["queryCoverage"])
        form.targetCoverage.data = Decimal(self.__data["targetCoverage"])
        form.eValue.data = Decimal(self.__data["eValue"])


    def save(
        self
    ):
        """
        Detailed description.
        """
        with open(self.PATH,"w") as ofile:
            ofile.write(dumps(self.__data))


    def update(
        self
        ,form
    ):
        """
        Detailed description.

        Parameters
        ----------
        form : 
        """
        self.__data = {
            "dataType": form.dataType.data
            ,"fpkm": str(form.fpkm.data)
            ,"singleEnd": form.singleEnd.data
            ,"seqComplete": form.seqComplete.data
            ,"minProteinLength": form.minProteinLength.data
            ,"refineStarts": form.refineStarts.data
            ,"outputFormat": form.outputFormat.data
            ,"taxonomy": form.taxonomy.data
            ,"queryCoverage": str(form.queryCoverage.data)
            ,"targetCoverage": str(form.targetCoverage.data)
            ,"eValue": str(form.eValue.data)
        }
