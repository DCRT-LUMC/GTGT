Usage
=====

.. _installation:

Installation
------------

gtgt can be installed using pip:

.. code-block:: console

   (.venv) $ pip install gtgt

To install the requirements for `gtgt server`, you can install from pip using:

.. code-block:: console

   (.venv) $ pip install gtgt[server]

Transcript
----------
You can fetch transcript information using the command below, or use the `/transcript` endpoint

.. code-block:: console

   gtgt transcript ENST00000241453.12 | jq .

Links
-----
You can fetch links to external resources for a specified variant using the command below, or use the `/links` endpoint

.. code-block:: console

   gtgt links "NM_002520.7:c.860_863dup"


Python functions
----------------
To work with Bed files, gtgt comes with a Bed class that will expand to BED12.

See :py:func:`gtgt.Bed`.
