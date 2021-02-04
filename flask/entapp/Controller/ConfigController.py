"""
Contains the ConfigController class.
"""




class ConfigController():
    """
    Detailed description.
    """
    PATH = "/workspace/entap/config.ini"


    def update(
        self
    ):
        """
        Detailed description.
        """
        if not pathExists(self.PATH):
            pass
        else:
            # check timestamps
