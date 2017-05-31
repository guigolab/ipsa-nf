FROM debian:stable

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "IPSA"

# install needed tools
RUN apt-get update --fix-missing -qq && apt-get install -y -q \
    curl \ 
    wget \
    gawk \
    build-essential \
    zlib1g-dev \
    libncurses5-dev  \
    libncursesw5-dev \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install maptools
RUN curl -fsSL https://github.com/pervouchine/maptools/archive/v3.2.tar.gz | tar xz && \
    mkdir -p maptools-3.2/bin && \
    make -C maptools-3.2 all

# install sjcount
RUN curl -fsSL https://github.com/pervouchine/sjcount/archive/v3.1.tar.gz | tar xz && \
    make -C sjcount-3.1 all

ENV PATH $PATH:/maptools-3.2/bin:/sjcount-3.1:/sjcount-3.1/samtools-0.1.18