Genetic Therapy Generator Toolkit
=================================

Introduction
------------
GTGT provides methods to assess whether a given pathogenic variant can be
corrected using genetic therapies (such as exon skipping). The unique feature
of GTGT is that it takes the actual sequence of the patient transcript into
account. In the case of exon skipping, rather than checking if the exon is
in-frame, GTGT will assess each exon to determine if skipping it will result in an
in-frame transcript in the patient. In our examples section, we have an
instance of a deletion in CLN3 where skipping an exon several exons upstream of
the deletion restores the reading frame and protein function.

Scope
-----
Determining whether a given patient can benefit from a genetic therapy is
extremely complex, from the cause and progression of the disease down to the
feasibility of a therapy. GTGT was intentionally designed to be narrow in scope,
so as to accurately answer a small part of this question: For a given variant,
which genetic therapy will restore most of the features of the wildtype
transcript.

Rather than taking a rule based approach, GTGT will simulate the effect of all
possible therapies. In the case of exon skipping, GTGT will simulate the effect
of skipping every exon except the first and the last (since these cannot
physically be skipped), and prioritize the therapies which are most likely to
benefit the patient. This approach allows GTGT to find therapeutic options
which are missed by other, rule based approaches to assessing variant
eligibility

Limitations
-----------
To assess the difference between the patient transcript and the transcript
modified by a therapy, GTGT uses a very straightforward measure. It simply
determines how many of the features lost in the patient are recovered in a
proposed therapy. The drawback of this approach is that it cannot take into
account the effect of the change, only its magnitude. This is most apparent
when GTGT is used to analyze a pathogenic missense mutation. GTGT counts this
change as the loss of only a single amino acid from the transcript. Any therapy
which skips this exon will lack much more than a single amino acid, so it will
also receive a lower score. Here, we rely on the user to interpret the proposed
therapies and determine which is the most suitable.

Splice altering mutations and/or intronic mutations cannot be interpreted by
GTGT, since there is no way to predict the effect of these kind of mutations on
the transcript. If you know the effect of your mutation, e.g. alternative
splicing, intron inclusion etc, you should provide these mutations to GTGT as
the input, instead of (only) the variant which causes the splicing defect.

Check out the :ref:`usage <usage>` section on how to use the website, or
:ref:`install <installation>` the project locally.

.. note::

   This project is under active development. GTGT has its documentation hosted on Read the Docs.

Contents
--------

.. toctree::

   usage
   transcript
   install
   api
   changelog
