Installation
============

.. _installation:

gtgt can be installed using pip:

.. code-block:: console

   (.venv) $ pip install gtgt

To install the requirements for `gtgt webserver`, you can install from pip using:

.. code-block:: console

   (.venv) $ pip install gtgt[webserver]

Configuration
-------------
Several configuration options can be set in GTGT using environmental variables. GTGT will check the following environmental variables:

.. list-table:: Environmental variables
  :header-rows: 1

  * - Variable
    - Description
  * - `GTGT_CACHE`
    - Full path to a folder to store the GTGT annotation cache
  * - `GTGT_VulExMap`
    - Full path to the data for vulnerable exons. Please contact the author to
      request a copy of this file for local usage
  * - `MUTALYZER_SETTINGS`
    - Full path to the Mutalyzer configuration file, which is used when retrieving transcript information


Mutalyzer configuration
~~~~~~~~~~~~~~~~~~~~~~~
The Mutalyzer configuration file is used when retrieving transcripts from NCBI and Ensembl. If you do not set this configuration, you are very likely to encounter errors and/or bans when using their API.

.. code-block:: text

  MUTALYZER_CACHE_DIR = /path/to/mutalyzer/cache
  MUTALYZER_FILE_CACHE_ADD = true
  EMAIL = <your_email_here>
  NCBI_API_KEY = <your_NCBI_API_KEY_here>

Transcript
----------
You can fetch transcript information using the command below.

.. code-block:: console

   gtgt transcript ENST00000241453.12

Links
-----
You can fetch links to external resources for a specified variant using the command below.

.. code-block:: console

   gtgt links "NM_002520.7:c.860_863dup"


Python functions
----------------
To work with Bed files, gtgt comes with a Bed class that will expand to BED12.

See :py:func:`gtgt.Bed`.
