=========
Changelog
=========

.. Newest changes should be on top.

.. This document is user facing. Please word the changes in such a way
.. that users understand how the changes affect the new version.

-------
v0.2.12
-------
+ Fix a bug which caused protein domains to not be updated

-------
v0.2.11
-------
+ Hotfix to roll back protein domains in the frontend

-------
v0.2.10
-------
+ Enable support for protein domains in the frontend

------
v0.2.9
------
+ Fix a bug with analyzing non-human Ensembl transcripts
+ Remove the API endpoints
+ Remove pydantic models which were only used by the API endpoints
+ Refactor Provider class (used to cache API calls)
+ Add support for querying UCSC tracks for a given Description


------
v0.2.8
------
+ Add support for refseq (NM) transcripts

------
v0.2.7
------
+ Fix a bug where the ExonViz figure cannot be drawn at the default scale
+ Drop support for python 3.9
+ Add support for python 3.14

------
v0.2.6
------
+ Fix circular import with exonviz

------
v0.2.5
------
+ Group Therapies together in the html
+ Show a figure of the input variants in the html
+ Add RNA and protein variant descriptions

------
v0.2.4
------
+ Add spinner to show progress to html
+ Add clickable examples to the html
+ Add ``gtgt export`` to the command line interface
+ Add separate mRNA and protein features to the Transcript
+ Extend documentation sections of the html

------
v0.2.3
------
+ Enable skipping two adjacent exons by default
+ Explicitly check that the variants are specified in ``c.`` format
+ Refactor internal representation of variants
+ Reduce analysis time

------
v0.2.2
------
+ Update the analysis results to give more information on proposed Therapies
+ Improve type annotations
+ Reduce the analysis time by 50% per exon
+ Fix a bug where the number of changed amino acids was overestimated

------
v0.2.1
------
+ Re-use mutalyzer Description object to reduce load on Ensembl API
+ Drop support for Python 3.8
+ Add support for Python 3.13
+ Add number of remaining basepairs to the transcript comparison
+ Rename default annotations for clarity
+ Fix a bug when the specified variant partially overlaps an exon
+ Fix a bug when VariantValidator is not available
+ Check user input for errors in web interface
+ Rename module from GTGT to gtgt

------
v0.2.0
------
+ Rename python package from GTGT to gtgt

------
v0.1.1
------
+ Add support for re-using mutalyzer Descriptions
+ Add support for small deletions
+ Fix a bug parsing VariantValidator payload
+ Add web interface using Flask

------
v0.1.0
------
+ Change the transcript endpoint to use POST and a TranscriptId model
+ Use genomic location for ClinVar link
+ Add support for Ensembl transcripts on the links endpoint
+ Refactor and simplify the BedModel
+ Validate HGVS is valid before querying Variant Validator
+ Change the links endpoint to use POST and HGVS
+ Add FastAPI endpoint for comparing transcripts
+ Add FastAPI endpoint for skipping exons
+ Add Flask web app
+ Drop support for Python 3.8
+ Add support for python 3.13

------
v0.0.2
------
+ Add FastAPI endpoint for transcripts and links
+ Add link to stringDB
+ Add link to gnomAD
+ Add API endpoint /links/ to fetch URL's for external references
+ Fetch knownGene track from UCSC
+ Fetch transcript information from Ensembl

------
v0.0.1
------
+ Add class for Transcripts
+ Validate that the input is valid Bed format
+ Add class for Bed format
+ Initial commit
