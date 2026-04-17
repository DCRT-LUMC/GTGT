Transcripts
===========

.. _transcript:

GTGT supports all RefSeq and Ensembl transcripts, from any organism. However,
not all features are available with every transcript, and protein domains are
only available for HG38.

There is a fundamental difference between Ensembl and RefSeq transcripts which
is important to understand when using GTGT.

Ensembl transcripts
-------------------
Ensembl transcripts (which start
with ENS) are defined as regions on the reference genome. This means that
Ensembl transcripts always have the exact same sequence as the reference genome
on which they are defined, and Ensembl transcripts also contain introns.

RefSeq transcripts
------------------
In contrast, RefSeq transcripts (which start with NM) are defined as the sequence
of the mature transcript, independent of a reference genome. This means that
the sequence of NM transcripts can differ from the reference genome.
Further more, NM transcripts are defined as mature transcripts, which means they do
not have introns. It is therefore not possible to specify intronic mutations on
an NM transcript, as is explained on the `HGVS website
<https://hgvs-nomenclature.org/stable/background/refseq/#DNAc>`_:

  "a coding DNA reference sequence does not contain intron or 5' and 3' gene
  flanking sequences and **can not be used to describe variants in introns** [...].
  Intronic sequences [...] may only be used to
  describe a variant when a genomic reference sequence identifier is provided."

Because NM transcripts do not have a relationship to the reference genome,
protein domain annotations from UCSC are not available for NM transcripts in
GTGT.

RefSeq transcripts on the genome
--------------------------------
Although NM transcript are in theory independent from the reference genome, in
practice they are closely related. To account for this, NM transcripts also
have an alternative version which is defined on a reference genome, which we
will refer to as NC(NM). For example for DMD, the transcript ``NM_004006.3`` is
annotated on the reference genome as ``NC_000023.11(NM_004006.3)``. This should
be read as "On sequence  ``NC_000023.11``, the transcript with name ``NM_004006.3``".
Intronic variants should always be specified on the NC(NM), and NC(NM)
transcripts should be used with GTGT if you want to make use of protein domain
information from UCSC. To convert an NM to NC(NM), you can use mutalyzer.
However, it is important to remember that the sequence of an NM and NC(NM)
transcript can be different!

An example of a RefSeq transcript that has a different sequence than the
reference chromosome is ``NM_017680.6``. If you go to the `mutalyzer mapper
<https://mutalyzer.nl/mapper?description=NM_017680.6%3Ac.%3D&reference_id=NC_000009.12&selector_id=NM_017680.6&slice_to=transcript>`_
you can see that ``NM_017680.6`` without any variants has three differences
with the transcript definition on the chromosome.

If you go to the regular `mutalyzer
<https://mutalyzer.nl/normalizer/NM_017680.6:c.100del>`_ website with your NM
transcript variant, you will see the option to get the `CHROMOSOMAL DESCRIPTIONS`,
with a warning banner if there are differences between the sequence of the
transcript and the reference genome.

MANE Select transcripts
-----------------------
MANE select transcripts are a special set of transcripts that act as the
default transcript for their respective genes, which are also guaranteed to
be identical between Ensembl and RefSeq transcripts. Since Ensembl transcripts
are guaranteed to be identical to the reference genome, this means that MANE
transcripts are also identical to the reference genome (NM and NC(NM) have the
exact same sequence).
