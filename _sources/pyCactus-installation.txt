Installation
------------
.. admonition:: Prerequisites
   :class: warning

   - `NumPy <https://pypi.python.org/pypi/numpy>`_ (tested with version 1.9.2)
   - `f90nml <https://pypi.python.org/pypi/f90nml>`_ (tested with version 0.12)
   - `PyEVTK <https://pypi.python.org/pypi/PyEVTK>`_ (tested with version 1.0) 

pyCactus may be installed by navigating to the `pyCactus` directory and executing:

.. code-block:: bash

	$ python setup.py install

If setuptools is correctly configured, this will also download and install module dependencies (`PyEVTK`, `f90nml`).

The :doc:`/scripts` will be automatically added to the system path.

Alternately, to install pyCactus in development mode (to be able to edit and test source files without having to reinstall):

.. code-block:: bash

	$ python setup.py develop