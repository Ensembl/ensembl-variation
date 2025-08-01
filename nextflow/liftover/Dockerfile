###############################################################################
# BCFtools +liftover
#   - Ubuntu 22.04 base
#   - Downloads latest stable BCFtools (1.21 as of 2025-05)
#   - Compiles freeseek/score plugins (inc. liftover)
###############################################################################

# --------------------------- Build -------------------------------------------
FROM ubuntu:22.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive

# build-time dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential wget curl ca-certificates \
        autoconf automake libtool pkg-config \
        zlib1g-dev libbz2-dev liblzma-dev \
        libcurl4-openssl-dev libssl-dev \
        libopenblas0-openmp libsuitesparse-dev \
    && rm -rf /var/lib/apt/lists/*

# fetch and unpack the latest bcftools tag
ARG BCFTOOLS_VERSION=1.21
WORKDIR /opt
RUN wget -q https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
 && tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2

# replace stock plugins with the freeseek/score versions
RUN set -eux; \
    cd bcftools-${BCFTOOLS_VERSION}/plugins; \
    for f in score.c score.h munge.c liftover.c metal.c blup.c; do rm -f "$f" || true; done; \
    wget -q -P . \
      https://raw.githubusercontent.com/freeseek/score/master/liftover.c

# compile bcftools + plugins
RUN cd bcftools-${BCFTOOLS_VERSION} \
 && make -j$(nproc) \
 && make install PREFIX=/usr/local \
 && mkdir -p /usr/local/lib/bcftools/plugins \
 && cp plugins/*.so /usr/local/lib/bcftools/plugins/

# --------------------------- Runtime -----------------------------------------
FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

# runtime libs
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4 libbz2-1.0 liblzma5 \
        libopenblas0-openmp libcholmod3 \
    && rm -rf /var/lib/apt/lists/*

# copy bcftools and its plugins from builder
COPY --from=builder /usr/local /usr/local
ENV BCFTOOLS_PLUGINS=/usr/local/lib/bcftools/plugins

ENTRYPOINT ["bcftools"]
CMD ["--help"]