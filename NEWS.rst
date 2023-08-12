Release Notes
=============

Changes in v 0.5.7:
*******************
- modified reporting of nt_diff so that it reports the number of nt differences between the query and the reference, including any length differences.
- updated the code to use Bio.Align rather than Bio.pairwise2, which is deprecated.

Changes in v 0.5.6:
*******************
- compatibility update for receptor_utils 0.0.40

Changes in v 0.5.5:
*******************
- speed optimisation for dig_sequence
- allow find_alignments to be called without a -ref parameter

Changes in V 0.5.4:
*******************
- Minor fixes to handling of non-functional sequences

Changes in V 0.5.3:
*******************
- Fixes to annotation in reverse-sense: false positives were not being filtered correctly

Changes in V 0.5.2:
*******************
- Added dig_sequence command, which allows a sequence stored locally or in Genbank to be searched for a specific allele. The closest match will be annotated.

Changes in V 0.5.1:
*******************
- Refactored elements of the code to make it more modular and easier to maintain.

Version 0.5: April 2023
***********************

First public version.

