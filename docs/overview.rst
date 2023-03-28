.. _overview_label:

Overview
========

Scope and Features
******************

Digger identifies receptor germline sequences in a genomic sequence or assembly. It requires two inputs: the sequence or assembly, and a set of reference genes to use as a starting point. An existing
reference set for the species under study is an ideal starting point, but if one is not available, good results can be obtained by using the reference set for another species, preferably one that is
reasonably closely related. The package is intended to be straightforward to use, but some familiarity with command-line tools is expected.

Digger starts by searching the sequence for approximate matches to sequences in the reference set. Where these are found, the match is extended to the full length of the matched sequence. A window
at either end is then checked for the expected flanking sequences (e.g. leader, RSS). These are identified by means of position weight matrices (PWMs). The flanking sequences have a well-established 
'canonical' form that is conserved between species. Digger contains PWMs for human and rhesus macaque IG loci, and these can be used as a starting point for other species. However, as some variation
is observed between species, Digger also contains tools for deriving tailored PWMs for a species of interest. These can be obtained from an existing annotation of the locus, or from an initial
annotation conducted with the human or rhseus PWMs.

A brief summary of the tools contained in the package is given below. In most cases, annotation can be performed with the `digger` tool. This will call subsidiary tools as necessary. Please refer 
to the Examples section for more information.

.. _FeatureTable:

.. csv-table::
   :file: tools/tool_summary.tsv
   :delim: tab
   :header-rows: 1
   :widths: 25, 75
   
