.. _parse_imgt_annotations:

parse_imgt_annotations
======================

``parse_imgt_annotations`` downloads an IMGT annotation file or or uses a file already downloaded. It parses the file to provide a list of annotated features.
Optionally it will also store the file downloaded, and create a FASTA file containing the annotated assembly.
Please refer to :ref:`human_igh` for example usage.

As of writing, IMGT has implemented a 'bot checker' mecchanism which means that automated downloads of annotation files are blocked. Unless this is modified in 
the future, you will need to download the annotation file manually from IMGT, and provide it to this tool using the ``-input_file`` argument. To do this:

1. Browse to the URL `https://www.imgt.org/ligmdb/view.action?format=IMGT&id=<id>` where ``<id>`` is the IMGT identifier of the sequence you wish to download annotations 
   for. For example, for the rhesus macaque IGH sequence with IMGT ID IMGT000064, the URL is `https://www.imgt.org/ligmdb/view.action?format=IMGT&id=IMGT000064`.

2. Check that you are on the IMGT View (see the list of near the top of the page, and click IMGT if necessary).

3. Press the Download button under the list of views (if you do not see the download button, check step 2).

Now proceed to run the ``parse_imgt_annotations`` tool with the ``-input_file`` argument set to the path of the file you just downloaded.

.. argparse::
   :filename: ../src/digger/parse_imgt_annotations.py
   :func: get_parser
   :prog: parse_imgt_annotations

