
***************
Getting Started
***************

Prerequisites
-------------

Before installation you need to ensure that your machine is equiped with the following:

   * `python 2.7 <http://www.python.org>`_ or later
   * `numpy <http://www.numpy.org>`_ python package
   * C++ compiler equipped with `C++11 <http://en.wikipedia.org/wiki/C%2B%2B11>`_ standard
   * MPI Library, for example `Open MPI <http://www.open-mpi.org>`_ or `MPICH <http://www.mpich.org>`_



Python
######

Please note that in order to install mapp4py you should have write privilages over the python that you are using. If not you can follow these instruction:

Download the latest version of `python <http://www.python.org>`_ and unpack it, enter the directory of source code and configure the source code. You might need to specify the location of installation to a place where you have write privilages. For example:

.. code-block:: bash

   $ ./configure --prefix=$HOME/MyPython

if the configuration finishes with out any problem you can compile the code by executing


.. code-block:: bash

   $ make

and to finalize the installation execute

.. code-block:: bash

   $ make install


NumPy
#####

Easiest way to install `numpy <http://www.numpy.org>`_ is by executing

.. code-block:: bash

   $ $HOME/MyPython/bin/pip install numpy

Please note that pip does not come with python installation by default. If your python is not equipped with pip you can follow the instruction `here <https://pip.pypa.io/en/stable/installing/#installing-with-get-pip-py>`_

.. code-block:: bash

   $ wget https://bootstrap.pypa.io/get-pip.py
   $ $HOME/MyPython/bin/python get-pip.py
   $ rm get-pip.py
   $ $HOME/MyPython/bin/pip install numpy




Download
--------

MAPP can be downloaded directly from `github <http://github.com/sinamoeini/mapp4py>`_ or simply execute the following in terminal

.. code-block:: bash

   $ git clone http://github.com/sinamoeini/mapp4py


Installation
------------

There are two ways that mapp4py can be installed, outlined below. My personal favorite is using Makefile due to its capability of compilation in parallel.

setup.py
########

Couple of the samples of setup.py are provided in the source code, choose the one that is more appropriate and edit it. The variables that need modifying are: :code:`mpi_lib`, :code:`mpi_lib_dir`, and :code:`mpi_include_dir`. If your cluster//computer is equiped with mpi compilers you can deduce these variables by executing

.. code-block:: bash

   $ mpic++ -show

If your default compilers are not equipped with C++11, you need to modify :code:`os.environ["CC"]` and :code:`os.environ["CXX"]`. Once the setup.py is ready execute:

.. code-block:: bash

   $ $HOME/MyPython/bin/python setup.py install


Makefile
########

Couple of the samples of Makefile are provided in the source code, choose the one that is more appropriate and edit it according to your cluster settings. To compile execute

.. code-block:: bash

   $ make py

If you want to compile in parallel instead you can execute:

.. code-block:: bash

   $ make py -j

to finallize installation:

.. code-block:: bash

   $ make install

Execution
---------

Onc you are done with installation and you have created an input script, mapp4py can be executed like regular python scripts. For example if you want to use 4 processors execute:

.. code-block:: bash

   $ mpiexec -n 4 $HOME/MyPython/bin/python input.py
