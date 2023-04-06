.. _compare_annotations:

compare_annotations
===================

``compare_annotations`` compares annotations produced by digger with IMGT's annotations as summarised by parse_annotations. Three files are produced:
- a .jpg showing Venn diagrams of identified functional annotations
- a .txt file summarising the specific functional sequences that were only identified in one annotation as opposed to both
- a .csv file listing agreements and differences of all sequences annotated by either mmethod
Please refer to :ref:`human_igh` for example usage.

.. argparse::
   :filename: ../src/digger/compare_annotations.py
   :func: get_parser
   :prog: compare_annotations

