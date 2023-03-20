.. _introduction_label:

Introduction
============

Digger is a toolkit for the automatic annotation of V,D and J genes and associated features in B cell receptor IGH genomic loci, developed for use with `VDJbase <https://vdjbase.org>`_. 

Pipeline
********

* Identify putative gene locations by running the assembly through BLAST, against an existing seed reference set, possibly from another species. BLAST must be `installed locally <https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_ (just the executable, no NCBI databases are needed).
* Process each high-likelihood match, identifying leader and RSS
* Extract confirmed gene sequences and features
* Optionally, SEARCH-D `(Safonova and Pevzner) <https://genome.cshlp.org/content/early/2020/10/15/gr.259598.119>`_ can be used to identify additional D genes de novo: this will give improved sensitivity for the short D-gene sequences. Please use this `forked version <https://github.com/williamdlees/SEARCH-D>`_, which has  minor changes so that gene co-ordinates are returned.
* Output is a CSV file of co-ordinates and features. Gapped sequences are provided using the closest sequence from a provided reference set.

In most cases, annotation can be performed by calling ``digger.py``. This will call other scripts in the package to run each step of the pipelinbe. These are documented in case they are needed individually in special cases.

Usage examples can be found in .bat files in the subdirectories. These should run in either Linux or Windows. 
