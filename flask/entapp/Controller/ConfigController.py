"""
Contains the ConfigController class.
"""
from ..Model.ConfigModel import *
from ..Model.ContamsModel import *
from ..Model.UninformsModel import *
from os import makedirs
from os.path import dirname
from os.path import getmtime




class ConfigController():
    """
    Detailed description.
    """
    PATH = "/workspace/entap/config.ini"
    HEADER_LINES = [
        "entap-db-bin=/workspace/entap/outfiles/bin/entap_database.bin"
        ,"entap-db-sql=/workspace/entap/outfiles/databases/entap_database.db"
        ,"entap-graph=/usr/local/bin/entap_graphing.py"
        ,"rsem-calculate-expression=rsem-calculate-expression"
        ,"rsem-sam-validator=rsem-sam-validator"
        ,"rsem-prepare-reference=rsem-prepare-reference"
        ,"convert-sam-for-rsem=convert-sam-for-rsem"
        ,"frame-selection=2"
        ,"genemarkst-exe=gmst.pl"
        ,"transdecoder-long-exe=/usr/lib/TransDecoder-v5.3.0/TransDecoder.LongOrfs"
        ,"transdecoder-predict-exe=/usr/lib/TransDecoder-v5.3.0/TransDecoder.Predict"
        ,"ontology=0,"
        ,"level=1,"
        ,"eggnog-sql=/workspace/entap/outfiles/databases/eggnog.db"
        ,"eggnog-dmnd=/workspace/entap/outfiles/bin/eggnog_proteins.dmnd"
        ,"interproscan-exe=/usr/local/interproscan/interproscan.sh"
        ,"protein="
        ,"diamond-exe=diamond"
    ]


    def update(
        self
    ):
        """
        Detailed description.
        """
        if not pathExists(dirname(self.PATH)):
            makedirs(dirname(self.PATH))
        if not pathExists(self.PATH):
            self.__write_()
        else:
            t1 = getmtime(ConfigModel.PATH)
            t2 = getmtime(ContamsModel.PATH)
            t3 = getmtime(UninformsModel.PATH)
            target = getmtime(self.PATH)
            if t1 > target or t2 > target or t3 > target:
                self.__write_()


    def __write_(
        self
    ):
        """
        Detailed description.
        """
        configModel = ConfigModel()
        contamsModel = ContamsModel()
        uninformsModel = UninformsModel()
        configs = self.HEADER_LINES
        configs += configModel.configLines()
        configs += contamsModel.configLines()
        configs += uninformsModel.configLines()
        with open(self.PATH,"w") as ofile:
            ofile.write("\n".join(configs) + "\n")
