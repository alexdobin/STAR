FROM debian:stretch-slim

MAINTAINER dobin@cshl.edu

ENV STAR_VERSION 2.7.0a

ENV PACKAGES gcc g++ make wget zlib1g-dev

RUN set -ex

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    cd /home && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz && \
    tar xzf ${STAR_VERSION}.tar.gz && \
    cd STAR-${STAR_VERSION}/source && \
    make STARstatic && \
    mkdir /home/bin && \
    cp STAR /home/bin && \
    cd /home && \
    'rm' -rf STAR-${STAR_VERSION} && \
    apt-get --purge autoremove -y  ${PACKAGES}

ENV PATH /home/bin:${PATH}

