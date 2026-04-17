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
CLN3 is a gene which is involved in the lysosome and loss of CLN3 is associated
with severe condition called Batten disease. In patients with Batten involving
CLN3, the most common mutation is a deletion of exon 8 and 9
(``ENST00000636147.2:c.461_677del``). `Click here
<https://gtgt.rnatherapy.nl/ENST00000636147.2:c.461_677del?protein_domains=true&extended=false>`__
to analyze this mutation using GTGT.

As can be seen from Transcript visualization, deletion of exon 8 and 9 causes a
frameshift mutation in the transcript. The frameshift leads to the loss of the
Lysosomal targetting motif, which is essential for the localisation of CLN3 to
the lysosome.

As can be seen from the figure, the reading frame can be restored by skipping
either the exon before the deletion (exon 7) or the exon after the deletion
(exon 10). Note that skipping exon 10 does not fully preserve the lysosomal targeting motif.

However, GTGT has found another therapy which restores the reading frame of the
transcript, namely skipping exon 6. This solution is not intuitive at all, and
also difficult to see from the transript visualization of CLN3. The best way to
explain it is that skipping exon 6 causes a frameshift when translating from
exon 5 to exon 7, which is then cancled by the frameshift caused by the exon
8+9 deletion.

Interestingly, exon 6 is a vulnerable exon that is already partially skipped in
healthy individuals. Skipping exon 6 has been shown to restore CLN3 function
and localisation to the lysosome in patient cell lines (`Centa2020
<https://www.nature.com/articles/s41591-020-0986-1>`__). Note that they use a
different transcript of CLN3 which lacks the first noncoding exon, so the
deletion is exon 7+8, and the exon they skip is exon 5.
