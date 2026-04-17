Examples
==========

.. _examples:

This page lists a number of examples of variants that can be analyzed using GTGT, sometimes with surprising results



Frameshift mutation in SDHD
---------------------------
The first example concerns a frameshift mutation in exon 2 of SDHD:
``ENST00000375549.8:c.149_150dup``. `Click here
<https://gtgt.rnatherapy.nl/ENST00000375549.8:c.149_150dup?protein_domains=true&extended=false>`__
to view the results in GTGT.

This variant is located in exon 2, which is an in frame exon. As expected,
skipping exon 2 is one of the therapies that is suggested by GTGT.

GTGT also lists skipping both exon 2 and 3 together, but as can be seen from
the annotations this is not a viable therapy. (This therapy was included in the
results due to the fact that it directly modifies the frameshift variant of the
patient).

However, GTGT also suggests a third therapy which will restore the reading
frame of the transcript, namely skipping exon 3. In the wildtype transcript,
exon 3 is not in an frame exon, and thus should not be eligible to skip.
However, since the patient as a frame shift mutation in exon 2, skipping exon 3
restores the reading frame again. Whether or not this will be beneficial for
the patient will depend on the exact features of the protein which are
preserved and the role they play in the disease. For example, we can see from
the results that exon 2 contains most of the mitochondrion trans membrane
peptide, which is lost when exon 2 is skipped. GTGT shows us that most of this
domain can be preserved if we skip exon 3 instead.

Premature STOP in DMD
---------------------
This example highlights how GTGT supports multiple variants on a single allele,
and correctly predict the resulting protein function:

``NC_000023.11(NM_004006.3):c.2504_2505delinsAG``. `Click here <https://gtgt.rnatherapy.nl/NC_000023.11(NM_004006.3):c.2504_2505delinsAG?protein_domains=true&extended=false>`__
to view the results in GTGT.

This example shows two mutations in exon 20 of DMD, which together introduce a
STOP codon. As can be seen from the transcript visualization in GTGT from
ExonViz, exon 20 is not an in frame exon, so it cannot be skipped to exclude
the introduced STOP codon.

However, we can also see from the results of GTGT that the combination of exon
19 and 20 and the combination of 20 and 21 are in frame, and lead to slightly
different losses of the Spectrin repeat 4 and 5. Determining which of these
therapies is most likely to benefit the patients in left to the user.


Multi-exon deletion in CLN3
---------------------------
CLN3 is a gene which is involved in the lysosome,
and the deletion of exon 8 and 9 is one of the most common causes for Batten
disease.

``ENST00000636147.2:c.461_677del``
