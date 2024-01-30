.. _docker:

Docker Image
----------------

The Docker image contains a working installation of Digger and its dependencies. The examples described in this documentation
are installed in the image at /digger/tests and can conveniently be run from the container.

The image is available in Docker Hub at `williamlees/digger <https://hub.docker.com/r/williamlees/digger>`_.

.. code-block:: bash

    # Pull the latest version
    $ docker pull williamlees/digger:latest
    
To use Digger from the command line, please use the following command in Linux/Mac/Windows:

.. code-block:: bash

    $ docker run williamlees/digger:latest digger --help
	

To run commands within the container, please use the following command in Linux/Mac:

.. code-block:: bash

    $ docker run -it -v $(pwd):/scratch williamlees/digger:latest bash

Or at the Windows command prompt:
	
.. code-block:: bat

    > docker run -it -v %cd%:/scratch  williamlees/digger:latest bash

This will start a bash shell in the container, with the current directory mounted at /scratch. You can then run Digger commands as described in the documentation.

