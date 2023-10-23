Usage
=====

.. _installation:

Installation
------------

GTGT can be installed using pip:

.. code-block:: console

   (.venv) $ pip install GTGT

Default usage
----------------
By default, Python-project will just print a greeting

.. code-block:: console

   python-project
   Hello, world

You can also specify the person to greeting

.. code-block:: console

   python-project --name John
   Hello, John

Python functions
----------------
To work with Bed files, GTGT comes with a Bed dataclass that will expand to BED12.

See :py:func:`GTGT.Bed`.
