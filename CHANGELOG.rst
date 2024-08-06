Changelog
=========

.. Newest changes should be on top.

.. This document is user facing. Please word the changes in such a way
.. that users understand how the changes affect the new version.

v0.0.3-dev
----------
+ Add support for Ensembl transcripts on the links endpoint
+ Refactor and simplify the BedModel
+ Validate HGVS is valid before querying Variant Validator
+ Change the links endpoint to use POST and HGVS
+ Add FastAPI endpoint for comparing transcripts
+ Add FastAPI endpoint for skipping exons

v0.0.2
------
+ Add FastAPI endpoint for transcripts and links
+ Add link to stringDB
+ Add link to gnomAD
+ Add API endpoint /links/ to fetch URL's for external references
+ Fetch knownGene track from UCSC
+ Fetch transcript information from Ensembl

v0.0.1
------
+ Add class for Transcripts
+ Validate that the input is valid Bed format
+ Add class for Bed format
+ Initial commit
