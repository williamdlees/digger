Release Notes
=============

Changes in 0.7.0
****************
- modify motif handling to allow multiple motifs in those cases where the length can be variable
- add variable length motifs for l-part1 and l-part2
- motifs are now based on the analysis in tests/utr
- add identification of TATA_BOX and OCTAMER in promoter regions

Changes in 0.6.9
****************
- bump dependency on receptor-utils to fix an issue with naming of D novel alleles

Changes in 0.6.8
****************
- improve documentation, add additional examples in tests, add Dockerised version
- include TRG motif files in the package

Changes in 0.6.7
****************
- include motif files for TRG

Changes in 0.6.7
****************
- better error handling and explanation in parse_imgt_annotations
- Docker image added

Changes in 0.6.6
****************
- fix crash in dig_sequence if annotating V-gene without specifying a gapped reference file
- update the package to include missing motif files

Changes in 0.6.4
****************
- fix issue in dig_sequence

Changes in 0.6.3:
*****************
- add support for multiple definitions of the J motif for a locus
- move J-TRP and other motif definitions to a file in the motifs directory

Changes in 0.6.2:
*****************
- require the user to supply an email address for requests to GenBank.

Changes in 0.6.1:
*****************
- bump required version of receptor-utils

Changes in 0.6.0:
*****************
- added support for TRA, TRB, TRD, TRG and test cases and motifs for human TR loci

Changes in v 0.5.10:
********************
- fixed unintended breakpoint in find_alignments which could fire under some circumstances

Changes in v 0.5.9:
*******************
- fixed erroneous reporting of leader1 sequence and co-ordinates when leader1 was not included in the query sequence

Changes in v 0.5.8:
*******************
- fixed various issues with end effects leading to the reporting of negative coordinates
- fixed an issue with targeted annotation which caused negative sense not to be properly reported
- changed IMGT urls to use https as http is no longer utilised

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

