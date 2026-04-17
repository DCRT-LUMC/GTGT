Using GTGT
==========

.. _usage:

GTGT is accessible online via `https://gtgt.rnatherapy.nl
<https://gtgt.rnatherapy.nl>`_. Users can specify their transcript and mutation
of interest in `HGVS <https://hgvs-nomenclature.org/stable/>`_ format. For the
best results, users should specify all variants that were observed on the
affected allele, since the exact patient sequence is taken into account when
assessing therapies. However, if only the pathogenic variant is known, this
will also work.

  .. figure:: figures/gtgt-input.png

     Input screen for the GTGT website. Users can describe their variant(s) of itnerest, and configure the parameters of the analysis.

Please refer to the :ref:`transcripts <transcript>` section for more
information on which transcripts can be used with GTGT.

Settings
--------
There are a number of settings that the user can change to control the analysis by GTGT.

Include protein domains
~~~~~~~~~~~~~~~~~~~~~~~
For transcripts encoded on HG38, GTGT can also access the protein features
which are available from the UCSC, and determine the effect of the therapy on
protein features. For example, GTGT will correctly determine if a protein
domain is lost due to an exon skip, and will also update the protein domain if
(part of) the domain has been lost due to a frameshift.

Show all evaluated therapies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GTGT will always evaluate all possible therapies to simulate their effect on the transcript. Select this option if you want to see output for all considered therapies.

By default, GTGT will only display a therapy when one of two conditions is met:

#. A therapy changes or removes at least one of the input variants
#. A therapy improves at least one of the features of the transcript, when compared to the input transcript

Option 2 can sometimes give surprising results, especially when protein domains are included in the analysis. If a therapy introduces a frameshift in the transcript, the region of the protein after the frameshift but before the STOP will sometimes contain motifs which match the wildtype protein. GTGT will interpret this as an improvement of that feature, and hence show this therapy in the results. In these cases, we rely on the users to recognize that the benefits of this therapy are extremely minimal.
