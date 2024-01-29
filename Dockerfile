FROM python:3.11-bullseye

MAINTAINER William Lees william@lees.org.uk

# Install BLAST+ executables
RUN BLAST=2.15.0 \
    && wget -q --show-progress --no-check-certificate \
       ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz \
    && tar -zxf ncbi-blast-${BLAST}+-x64-linux.tar.gz \
    && mv ncbi-blast-${BLAST}+/bin/* /usr/local/bin \
    && rm -r ncbi-blast-${BLAST}+-x64-linux.tar.gz ncbi-blast-${BLAST}+

# Install Digger
RUN \
  apt-get update && \
  apt-get -y --no-install-recommends install apt-utils && \
  apt-get -y install vim && \
  apt-get -y install git && \  
  git clone "https://github.com/williamdlees/digger" /digger && \
  pip install receptor-digger
  
