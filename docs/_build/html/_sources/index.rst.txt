.. Digger documentation master file

Digger: tools for annotating genomic assemblies of the IG/TR receptor loci
============================================================================
Digger is a toolkit for the automatic annotation of unrearranged V,D and J genes in B- and T- cell immunoglobulin receptor genomic loci. It can be used to annotate both entire assemblies, large fragments of a locus, 
or small fragments. It annotates all features of the gene (e.g. leader, RSS) excluding UTR, if the features are present. It attempts to classify the gene as Functional,
ORF, or pseudogene following IMGT practce.


.. toctree::
  :maxdepth: 1
  :caption: Getting Started

  overview
  docker	
  install
  news
  
  
.. toctree::
  :maxdepth: 1
  :caption: Examples

  examples/human_igh
  examples/rhesus_igh
  examples/targeted_annotation
  examples/additional_examples
  
.. toctree::
  :maxdepth: 1
  :caption: Usage Documentation

  usage
  tools/annotation
  

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
