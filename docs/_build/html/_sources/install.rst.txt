.. _install:



Installation
============

The latest stable version of Digger may be downloaded from `PyPI <https://pypi.python.org/pypi/receptor-digger>`__. Digger requires Python 3.9 or above.

Development versions are available from `GitHub <https://github.com/williamdlees/digger>`__.

Digger requires `BLAST <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`__ to be installed. If you use conda/anaconda on Linux or Mac we recommend installing BLAST with conda:

    > conda install bioconda::blast

Alternatively, please follow the link for installation instructions. No BLAST databases are needed: just the executable. 

Once BLAST has been installed, please verify by typing `blastn --help` at the command line: if everything is ok, it should provide
usage instructions.

If you encounter the error ``Cannot allocate memory``, this is coming from BLAST makeblastdb. Please set the environment variable ``BLASTDB_LMDB_MAP_SIZE=100000000``. 

If you encounter the error ``Error while loading shared libraries: libnsl.so.1``, the shared library ``libnsl.so.1`` needs to be installed. The precise instructions to do so are platform specific, but you should be able to obtain guidance through a Google search. 


The easiest way to install Digger itself is with pip::

    > pip install receptor-digger --user

Digger requires the Python packages `receptor-utils <https://pypi.org/project/receptor-utils/>`__ and `biopython <https://pypi.org/project/biopython/>`__. These should be installed automatically by pip.

