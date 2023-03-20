.. _digger_label:

digger
======

``digger`` annotates a single reference assembly, using BLAST to search the assembly for potential germline sequences. It requires an initial reference set
for BLAST to use: this could come from a similar species, or a former annotation.

.. argparse::
   :filename: ../src/digger/digger.py
   :func: get_parser
   :prog: digger


At least one file containing reference genes must be provided. You can, for example, supply ``v_ref``, ``d_ref`` and ``j_ref``, or just ``v_ref``. Digger will annotate whatever genes are discovered with the corresponding set(s). 
In practice, the sets do not have to be that good a match: BLAST will identify partial matches, and Digger's logic will extend the match to a full gene, including canonical RSS and leader (using the ``motif`` folder).

``v_ref_gapped`` is used to gep v-sequences correctly in order to identify conserved codons and so on. Again these sequences do not need to be that good a match in practice. The sequences **must be IMGT aligned with
no extraneous codons**. Note in particular that IMGT has introduced insertions into macaque alignments in recent years. **Sets with these insertions should not be used**.

``ref_comp`` allows you to specify that you would like annotated sequences to be compared with sequences in a set. You can include as many different sets as you wish. The output file will contain columns
for each of these, listing the closest sequence found and the proximity (%, and number of nucleotides).

If you choose not to specify the ``sense``, Digger will select the sense that elicits the highest evalue (results are shown in the output so that you can decide whether it has made the right choice, 
and whether you wish to annotate in both senses)

If ``searchd`` is specified, Digger will infer D genes using SEARCH-D rather than BLAST. In this case, ``d_ref`` is not used. ``searchd`` **must specify the path to the forked version of SEARCHD** at x. This generates
co-ordinates for each gene, which are not provided by the original code.
