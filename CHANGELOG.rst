Changelog
=========

.. Newest changes should be on top.

.. This document is user facing. Please word the changes in such a way
.. that users understand how the changes affect the new version.

v0.2.4
------
+ Add spinner to show progress to html
+ Add clickable examples to the html
+ Add ``gtgt export`` to the command line interface
+ Add separate mRNA and protein features to the Transcript
+ Extend documentation sections of the html

v0.2.3
------
+ Enable skipping two adjacent exons by default
+ Explicitly check that the variants are specified in ``c.`` format
+ Refactor internal representation of variants
+ Reduce analysis time

v0.2.2
------
+ Update the analysis results to give more information on proposed Therapies
+ Improve type annotations
+ Reduce the analysis time by 50% per exon
+ Fix a bug where the number of changed amino acids was overestimated

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

v0.2.0
------
+ Rename python package from GTGT to gtgt

v0.1.1
----------
+ Add support for re-using mutalyzer Descriptions
+ Add support for small deletions
+ Fix a bug parsing VariantValidator payload
+ Add web interface using Flask

v0.1.0
----------
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
