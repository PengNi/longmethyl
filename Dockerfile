# Set the base image to Ubuntu and NVIDIA GPU from https://hub.docker.com/r/nvidia/cuda
# or from https://ngc.nvidia.com/catalog/containers/nvidia:cuda/tags

# CUDA 10.0.130 >= 410.48 (linux)  >= 411.31 (windows)
FROM nvidia/cuda:10.0-cudnn7-runtime-ubuntu18.04
# FROM nvidia/cuda:10.0-cudnn7-devel-ubuntu18.04

# Author and maintainer
MAINTAINER Peng Ni <543943952@qq.com>
LABEL description="longmethyl" \
      author="543943952@qq.com"

ARG DNAME="longmethyl"

# shouldn't do this?
# ENV CUDA_HOME /usr/local/cuda
# ENV PATH ${CUDA_HOME}/bin:$PATH
# ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:${CUDA_HOME}/compat/:${CUDA_HOME}/lib64/

# Guppy version
ARG GUPPY_VERSION=4.2.2  # 4.2.2, lastest version for cuda-10. However, maybe the guppy version is 
                         # not associated with the cuda in docker - it is associated with the 
                         # cuda/driver in the host machine?
ARG BUILD_PACKAGES="wget apt-transport-https procps git curl"
ARG DEBIAN_FRONTEND="noninteractive"

# Install guppy-gpu version, ref: https://github.com/GenomicParisCentre/dockerfiles
# trial 1 (docker: cpu-T, cpu-on-gpu-T, gpu-on-gpu-F; singularity: cpu-T, cpu-on-gpu-T, gpu-on-gpu-F)
RUN apt-get -q update && \
    DEBIAN_FRONTEND=${DEBIAN_FRONTEND} apt-get -q install --yes ${BUILD_PACKAGES} \
    libnvidia-compute-418 \
    libnvidia-common-418 \
    nvidia-dkms-418 && \
    cd /tmp && \
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb && \
    DEBIAN_FRONTEND=${DEBIAN_FRONTEND} apt-get -q install --yes /tmp/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb && \
    rm *.deb && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# trial 2 (docker: cpu-F, cpu-on-gpu-t, gpu-on-gpu-T; singularity: cpu-F, cpu-on-gpu-T, gpu-on-gpu-T)
# RUN apt-get -q update && \
#     DEBIAN_FRONTEND=${DEBIAN_FRONTEND} apt-get -q install --yes ${BUILD_PACKAGES} \
#           libcurl4 \
#           libssl-dev \
#           libhdf5-cpp-100 \
#           libzmq5 \
#           libboost-atomic1.65.1 \
#           libboost-chrono1.65.1 \
#           libboost-date-time1.65.1 \
#           libboost-filesystem1.65.1 \
#           libboost-program-options1.65.1 \
#           libboost-regex1.65.1 \
#           libboost-system1.65.1 \
#           libboost-log1.65.1 \
#           libboost-iostreams1.65.1 && \
#     cd /tmp && \
#     wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb && \
#     # dpkg -i --ignore-depends=nvidia-384,libcuda1-384 /tmp/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb && \
#     dpkg -i /tmp/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb && \
#     rm *.deb && \
#     apt-get autoremove --purge --yes && \
#     apt-get clean && \
#     rm -rf /var/lib/apt/lists/*

# trial 3
# RUN apt-get -q update && \
#    DEBIAN_FRONTEND=${DEBIAN_FRONTEND} apt-get -q install --yes ${BUILD_PACKAGES} && \ 
#    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPY_VERSION}_linux64.tar.gz && \
#    tar zxvf /tmp/ont-guppy_${GUPPY_VERSION}_linux64.tar.gz -C /opt && \
#    rm *.tar.gz && \
#    apt-get autoremove --purge --yes && \
#    apt-get clean && \
#    rm -rf /var/lib/apt/lists/*
# ENV PATH /opt/ont-guppy/bin:$PATH

# set vbz plugin
ARG VBZ_VERSION=1.0.1
RUN cd /tmp && \
    wget -q https://github.com/nanoporetech/vbz_compression/releases/download/v${VBZ_VERSION}/ont-vbz-hdf-plugin-${VBZ_VERSION}-Linux-x86_64.tar.gz && \
    tar zxvf ont-vbz-hdf-plugin-${VBZ_VERSION}-Linux-x86_64.tar.gz -C /opt && \
    rm *tar.gz
ENV HDF5_PLUGIN_PATH /opt/ont-vbz-hdf-plugin-${VBZ_VERSION}-Linux/usr/local/hdf5/lib/plugin

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
