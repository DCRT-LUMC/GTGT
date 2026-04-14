Use the GTGT website
====================

.. _usage:

GTGT is accessible online via `https://gtgt.rnatherapy.nl <https://gtgt.rnatherapy.nl>`_.
Users can specify their transcript and mutation of interest in HGVS format. For
the best results, users should specify all variants that were observed on the
affected allele, since the exact patient sequence is taken into account when
assessing therapies. However, if only the pathogenic variant is known, this
will also work. For transcripts encoded on HG38, GTGT can also access the
protein features which are available from the UCSC, and determine the effect of
the therapy on protein features. For example, GTGT wil correctly determine if a
protein domain is lost due to an exon skip, and will also update the protein
domain if (part of) the domain is frameshifted.
