# Set the base image to Ubuntu 18.04 and NVIDIA GPU from https://hub.docker.com/r/nvidia/cuda
# or from https://ngc.nvidia.com/catalog/containers/nvidia:cuda/tags

# FROM nvidia/cuda:10.0-base-ubuntu18.04
FROM nvidia/cuda:10.0-cudnn7-runtime-ubuntu18.04

# Author and maintainer
MAINTAINER Peng Ni <543943952@qq.com>
LABEL description="longmethyl" \
      author="543943952@qq.com"

# Guppy version
ARG DNAME="longmethyl"
ARG GUPPY_VERSION=4.2.2
ARG BUILD_PACKAGES="wget apt-transport-https procps git curl libnvidia-compute-450-server"
ARG DEBIAN_FRONTEND="noninteractive"

# Install guppy-gpu version, ref: https://github.com/GenomicParisCentre/dockerfiles
RUN apt-get -q update && \
    DEBIAN_FRONTEND=${DEBIAN_FRONTEND} apt-get -q install --yes ${BUILD_PACKAGES} && \
    cd /tmp && \
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb && \
    DEBIAN_FRONTEND=${DEBIAN_FRONTEND} apt-get -q install --yes /tmp/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb && \
    rm *.deb && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

#Install miniconda
RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda.sh && \
    /bin/bash Miniconda.sh -b -p /opt/conda && \
    rm Miniconda.sh

# Adding conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Create the environment:
COPY environment.yml /
RUN conda env create --name ${DNAME} --file=environment.yml && conda clean -a

# Make RUN commands use the new environment
# naem need to be the same with the above ${DNAME}
SHELL ["conda", "run", "-n", "longmethyl", "/bin/bash", "-c"]

# Install latest version for megalodon, even conflicts with fast5mod, they can work
# RUN pip install megalodon==${MEGALODON_VERSION} &&\
# 	pip install ont-remora==${REMORA_VERSION} &&\
#     pip cache purge &&\
#     npm install -g inliner && npm cache clean --force

# Set env path into PATH
ENV PATH /opt/conda/envs/${DNAME}/bin:$PATH
USER root
WORKDIR /data/

# Get bigwig conversion tool
# RUN cd /data && \
#     wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig &&\
#     chmod +x bedGraphToBigWig &&\
#     mv bedGraphToBigWig  /usr/local/bin/

RUN cd /data

CMD ["bash"]
