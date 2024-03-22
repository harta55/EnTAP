Installing EnTAP
===============================

EnTAP is packaged with all of the software necessary to fully annotate a set of transcripts.  It is optimized to allow a single-command execution for all steps in the pathway, including paramterization by the user.  EnTAP does not have a graphical user interface but it does generate visual summaries for the user at each stage as well as detailed summary files and logs. EnTAP must be installed and configured in order to begin annotating! A test dataset comes with EnTAP to ensure it has been configured properly.

EnTAP may be installed through the source code or the dockerfile. After installation is complete, Configuration must be ran once in order to download any necessary databases. Before moving to Configuration, please review the ini files first!

.. toctree::
   :maxdepth: 3

   Installation_From_Source_Code/installation_from_source_code.rst
   Installation_From_Dockerfile/installation_from_dockerfile.rst
