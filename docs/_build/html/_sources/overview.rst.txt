.. _overview_label:

Overview
========

There are two principal tools: ``digger`` itself, and ``dig-sequence``. ``digger`` performs de-novo annotation, searching for all germline sequences that are similar
to those in a reference set. ``dig_sequence`` searches for the single sequence that is most similar to a given reference sequence. Both tools use the same 
underlying code. ``dig_sequence`` can be useful for determining the full annotation of a genomic sequence you have identified, or for verifying a sequence
discovered in GenBank. ``digger`` is slower, but performs the full annotation of a genomic sequence or assembly.

digger
======

Scope and Features
******************

``digger`` identifies immunoglobulin receptor germline sequences in a genomic sequence or assembly. It requires two inputs: the sequence or assembly, and a set of reference genes to use as a starting point. An existing
reference set for the species under study is an ideal starting point, but if one is not available, good results can be obtained by using the reference set for another species, preferably one that is
reasonably closely related. The package is intended to be straightforward to use, but some familiarity with command-line tools is expected.

``digger`` starts by searching the sequence for approximate matches to sequences in the reference set. Where these are found, the match is extended to the full length of the matched sequence. A window
at either end is then checked for the expected flanking sequences (e.g. leader, RSS). These are identified by means of position weight matrices (PWMs). The flanking sequences have a well-established 
'canonical' form that is conserved between species. ``digger`` contains PWMs for human IG/TR and rhesus macaque IG loci, and these can be used as a starting point for other species. However, as some variation
is observed between species, The toolkit also contains tools for deriving tailored PWMs for a species of interest. These can be obtained from an existing annotation of the locus, or from an initial
annotation conducted with the human or rhesus PWMs. Digger's built-in PWM's are based on analysis of selected IMGT annotations.

Pipeline
********

.. image:: pipeline.jpg
   :width: 600

A brief summary of the tools contained in the package is given below. In most cases, annotation can be performed with the :ref:`digger` tool. This will call subsidiary tools as necessary. Please refer 
to the Examples section for more information. The first example, :ref:`human_igh`, explains straightforward use. The second example, :ref:`rhesus_igh`, shows how
to generate PWMs for a locus or species that is not already covered by the package. It also shows how to call the underlying tools directly.

dig-sequence
============

``dig_sequence`` identifies and annotates the closest match to a single specified reference sequence. It provides three modes: one in which the sequence to search
is provided in a FASTA file, one in which it is specified as a GenBank accession number, and a third in which multiple GenBank accession numbers to search are provided 
in a file. Its use is described in :ref:`dig_sequence`.

.. _FeatureTable:

.. csv-table::
   :file: tools/tool_summary.tsv
   :delim: tab
   :header-rows: 1
   :widths: 25, 75
   
