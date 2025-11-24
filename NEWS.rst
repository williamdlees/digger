Release Notes
=============

CHanges in 0.8.1
================
- Revised motif tables to fully implement the new Max/Min values for RSS spacers

Changes in 0.8.0
================
- Breaking change: the motif file 'motif_params.json' now takes minimum and maximum RSS spacings for V, D and J RSS, to allow tuning/exploration. 
  If you have a custom motif_params.json file, you will need to update it to the new format. Note that, by default, Digger has always provided 1nt
  of flexibility in V and J RSS spacing, eg 22-23, 11-12. This was previously built into the code, but now needs to be specified in the motif_params.json file.
  For examples of the new format, see `github <https://github.com/williamdlees/digger/tree/main/src/digger/motifs>`_.
- Added new exon annotation columns: exon1, exon1_start, exon1_start_rev, exon1_end, exon1_end_rev and corresponding columns for exon 2
- Added new RSS annotation columns: v_spacer, v_spacer_length, j_spacer, j_spacer_length, d_3_spacer, d_3_spacer_length, d_5_spacer, d_5_spacer_length
- Added new leader splice columns: donor-splice, acceptor-splice
- Revised annotation of L-PART1 to match IMGT convention (see :ref:`leader_annotation` for details). 
- Breaking change: REMOVED coordinate columns for l-part1, l-part2 - see above for why.

Changes in 0.7.7
================
- The donor splice (GT) is no longer included in the annotation of L-PART1. This is a break in custom with IMGT, who typically include the GT in the L-PART1 annotation. 
  However, this means that L-PART1 and L-PART2 now exactly match the exon, which we feel is important for clarity. To enable the donor splice to be examined more easily,
  Digger now explicitly annotates the V intron.
- Fix erroneous 'failed coordinate check' messages which were emitted where an annotated region extends beyond the assembly boundaries.

Changes in 0.7.6
================
- fix problem in handling multiple contigs

Changes in 0.7.5
================
- further salmonid name fix

Changes in 0.7.4
================
- fix problem with salmonid names, e.g. TRB3V

Changes in 0.7.3
================
- fix crash in compare_annotations.py

Changes in 0.7.2
****************
- digger will default to using the starting set of reference genes for comparison, if no sets are specified with -ref_comp.

Changes in 0.7.1
****************
- digger will now copy reference set files to the work directory if necessary, as blast's makedb requires the files to be local.
- digger can now handle multiple sequences (e.g. contigs) in the assembly file.

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

