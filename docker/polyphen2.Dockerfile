# Linux/x86_64 is required -- emulated in ARM platforms (e.g. newer Mac models)
FROM --platform=linux/x86_64 bioperl/bioperl

ARG pph2_version=polyphen-2.2.3r407
ARG blast_version=ncbi-blast-2.12.0+

# Install dependencies
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y \
  wget \ 
  build-essential \
  # Java RE
  default-jre \
  # Perl dependencies
  libscalar-list-utils-perl \
  libxml-simple-perl \
  libdbd-sqlite3-perl \
  libcgi-pm-perl \
  && rm -rf /var/lib/apt/lists/*

# Download PolyPhen-2
WORKDIR /opt
ENV PPH_DIR /opt/pph2
ENV PPH_BIN $PPH_DIR/bin
ENV PPH_TAR $pph2_version.tar.gz
RUN wget http://genetics.bwh.harvard.edu/pph2/dokuwiki/_media/$PPH_TAR && \
    tar vxaf $PPH_TAR && \
    rm $PPH_TAR && \
    mv polyphen* $PPH_DIR

# Install BLAT
WORKDIR $PPH_BIN
#RUN apt-get install -y curl
ENV BLAT_URL https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64
RUN wget $BLAT_URL/blat/blat $BLAT_URL/twoBitToFa $BLAT_URL/bigWigToWig && \
    chmod +x *

# Install BLAST
WORKDIR $PPH_DIR
ENV BLAST_DIR $PPH_DIR/blast
ENV BLAST_TAR $blast_version-x64-linux.tar.gz
ENV BLAST_URL https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+
# Parse version of BLAST to get proper URL directory
RUN BLAST_VERSION=$( echo $blast_version | grep -o "[0-9].*[0-9]" ) && \
    wget "$BLAST_URL/$BLAST_VERSION/$BLAST_TAR" && \
    tar xzvf $BLAST_TAR && \
    rm $BLAST_TAR && \
    rmdir $BLAST_DIR && \
    mv $blast_version $BLAST_DIR && \
    chmod +x $BLAST_DIR/*
ENV PATH "$PATH:$BLAST_DIR/bin"

# Build PolyPhen-2
WORKDIR $PPH_DIR/src
RUN apt-get update && apt-get install -y unzip && \
    # Avoid checking expired SSL certificates
    echo check_certificate = off > ~/.wgetrc && \
    make download && \
    # Uninstall unzip and remove global wget file
    apt-get install unzip && rm -rf /var/lib/apt/lists/* && rm ~/.wgetrc && \
    make clean && \
    make && \
    make install

WORKDIR $PPH_DIR
RUN ./configure
ENV PATH "$PATH:$PPH_BIN"

# Move data directories to /opt/pph2/data and symlink as expected by PolyPhen-2
# Polyphen-2 data can thus be mounted with a single bind mount, e.g.:
#   docker run -v ${HOME}/pph2_data:/opt/pph2/data -ti pph2
RUN mkdir data && \
    mv dssp nrdb models pdb2fasta precomputed ucsc uniprot wwpdb data && \
    ln -s data/* .
