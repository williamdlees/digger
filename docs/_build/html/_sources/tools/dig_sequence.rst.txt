.. _dig_sequence:

dig_sequence
============

``dig_sequence`` provides targeted annotation of genomic sequences. The sequences can be as small as a single coding region,
or as large as an entire locus. The tool will search for the sequence that best matches a specified target sequence, and
annotate just that best match. The target sequence is specified as a sequence name in a FASTA file: typically this would
be a germline set.

Options are available to annotate a single sequence in a FASTA file, a single genbank ID (the sequence will be fetched from
Genbank), or a list of sequences or Genbank IDs specified in a CSV file.


.. argparse::
   :filename: ../src/digger/dig_sequence.py
   :func: get_parser
   :prog: dig_sequence

