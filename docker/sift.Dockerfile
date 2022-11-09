# Linux/x86_64 is required -- emulated in ARM platforms (e.g. newer Mac models)
FROM --platform=linux/x86_64 ubuntu:20.04

ARG sift_version=sift6.2.1
ARG blast_version=ncbi-blast-2.12.0+

RUN apt-get update && apt-get -y install \
  wget \
  # csh is required to run SIFT scripts
  csh \
  && rm -rf /var/lib/apt/lists/*

# Install BLAST
WORKDIR /opt
ENV BLAST_TAR $blast_version-x64-linux.tar.gz
ENV BLAST_URL https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+
ENV BLAST_TARGET_DIR /opt/blast
# Parse version of BLAST to get proper URL directory
RUN BLAST_VERSION=$( echo $blast_version | grep -o "[0-9].*[0-9]" ) && \
    wget "$BLAST_URL/$BLAST_VERSION/$BLAST_TAR" && \
    tar xzvf $BLAST_TAR && \
    rm $BLAST_TAR && \
    mv $blast_version $BLAST_TARGET_DIR && \
    chmod +x $BLAST_TARGET_DIR/*

# Install SIFT
WORKDIR /opt
ENV TAR $sift_version.tar.gz
ENV SIFT_DIR /opt/sift
RUN wget https://s3.amazonaws.com/sift-public/nsSNV/$TAR && \
    tar -xzvf $TAR && \
    rm $TAR && \
    mv sift* $SIFT_DIR

# Set environmental (hardcoded) paths in SIFT_for_submitting_*.csh files
WORKDIR $SIFT_DIR/bin
ENV BLAST_BIN $BLAST_TARGET_DIR/bin
ENV BLIMPS_DIR $SIFT_DIR/blimps
RUN sed -Ei "s|(setenv NCBI) .*|\1 $BLAST_BIN|g" SIFT_for_submitting_*.csh
RUN sed -Ei "s|(setenv SIFT_DIR) .*|\1 $SIFT_DIR|g" SIFT_for_submitting_*.csh
RUN sed -Ei "s|(setenv BLIMPS_DIR) .*|\1 $BLIMPS_DIR|g" SIFT_for_submitting_*.csh

# Improve seqs_chosen_via_median_info.csh
ENV seqs_chosen seqs_chosen_via_median_info.csh
RUN sed -Ei 's|\$NCBI/||g' $seqs_chosen
RUN sed -Ei 's|.*(set tmpdir =).*|\1 "."|g' $seqs_chosen
# RUN sed -Ei 's|(.*rm \$pid\.TEMP\*)|#\1|g' $seqs_chosen
RUN sed -Ei 's|(.*limit cputime.*)|#\1|g' $seqs_chosen # avoid CPU time limit
RUN sed -Ei 's|(set median_threshold = )2.75;|\1$3|g' $seqs_chosen
RUN sed -i '/exiting/{n;s/exit (-1)/#exit (-1)/}' $seqs_chosen

WORKDIR $SIFT_DIR
ENV PATH "$PATH:$SIFT_DIR/bin"
ENV PATH "$PATH:$BLAST_BIN"
