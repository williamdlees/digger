.. _digger:

digger
======

``digger`` annotates a single reference assembly, using BLAST to search the assembly for potential germline sequences. It requires an initial reference set
for BLAST to use: this could come from a similar species, or a former annotation.
Please refer to :ref:`human_igh` for example usage.

.. argparse::
   :filename: ../src/digger/digger.py
   :func: get_parser
   :prog: digger


At least one file containing reference genes must be provided. You can, for example, supply ``v_ref``, ``d_ref`` and ``j_ref``, or just ``v_ref``. Digger will annotate whatever genes are discovered with the corresponding set(s). 
In practice, the sets do not have to be that good a match: BLAST will identify partial matches, and Digger's logic will extend the match to a full gene, including canonical RSS and leader (using the ``motif`` folder).

Digger requires a set of postion-weighted matrices, to identify RSS and leader. It is also possible to specify conserved locations of motifs. This `motif` data should be stored in a ``motif`` folder. Motifs for
human and rhesus macaque IG are built in to the package, and may be used with ``-species`` by specifying either ``human`` or ``rhesus_macaque``. The species is used in conjunction with ``-locus`` to determine
the correct motifs. Alternatively, ``-motif_dir`` can be used to specify custom motifs created outside of the package. Please refer to :ref:`calc_motifs` and to :ref:`rhesus_igh` for further details
on custom motifs.

``v_ref_gapped`` is used to gep v-sequences correctly in order to identify conserved codons and so on. Again these sequences do not need to be that good a match in practice. The sequences **must be IMGT aligned with
no extraneous codons**. Note in particular that IMGT has introduced insertions into macaque alignments in recent years. **Sets with these insertions should not be used**.

``ref_comp`` allows you to specify that you would like annotated sequences to be compared with sequences in a set. You can include as many different sets as you wish. The output file will contain columns
for each of these, listing the closest sequence found and the proximity (%, and number of nucleotides).

If you choose not to specify the ``sense``, Digger will select the sense that elicits the highest number of hits and the highest evalue (results are shown in the output so that you can decide whether it has made the right choice, 
and whether you wish to annotate in both senses)

