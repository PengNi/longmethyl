# longmethyl

__Learning Nextflow__ - A demo nextflow pipeline of methylation detection using long reads

## Installation

  - (1) Install conda from [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) if neeeded.


  - (2) Install nextflow.

Create an environment containing nextflow/install nextflow:
```sh
# create a new environment and install nextflow in it
conda create -n nextflow -c conda-forge -c bioconda nextflow

# or install nextflow in an existing environment
conda install -c conda-forge -c bioconda nextflow
```

  - (3) Download longmethyl from github.

```sh
git clone https://github.com/PengNi/longmethyl.git
```

  - (4) [optional] Install graphviz.

```sh
sudo apt install graphviz
# or
sudo yum install graphviz
```

## Usage

### Option 1. Run with singularity (recommend)

If it is the first time you run with singularity (e.g. using `-profile singularity`), the following cmd will cache the dafault singularity image (`--singularity_name`) to the `--singularity_cache` directory (default: `local_singularity_cache`) first. There will be a `.img` file in the `--singularity_cache` directory.

```sh
# activate nextflow environment
conda activate nextflow

# run longmethyl, this cmd will cache a singularity image before processing the data
nextflow run ~/tools/longmethyl -profile singularity \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
# or, run longmethyl using GPU, set CUDA_VISIBLE_DEVICES
CUDA_VISIBLE_DEVICES=0 nextflow run ~/tools/longmethyl -profile singularity \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
```

The downloaded `.img` file can be used then, without downloaded again:

```sh
# this time nextflow will not download the singularity image again, it has already
# been in the --singularity_cache directory.
nextflow run ~/tools/longmethyl -profile singularity \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
# or
nextflow run ~/tools/longmethyl -profile singularity \
    --singularity_cache local_singularity_cache \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
# or
nextflow run ~/tools/longmethyl -profile singularity \
    --singularity_name local_singularity_cache/nipengcsu-longmethyl-0.1.img \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
```

The singularity image can be also pulled before running the cmd. The pulled `.sif` file is only needed to be downloaded once.

```sh
# pull singularity image (once for all). There will be a .sif file. 
singularity pull docker://nipengcsu/longmethyl:0.1

# run longmethyl
nextflow run ~/tools/longmethyl -profile singularity \
    --singularity_name longmethyl_0.1.sif \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
```

### Option 2. Run with docker

  - (1) Pull docker image (once for all).

It is better to pull docker image before running pipeline the first time, cause this may be time-consuming and there may be network issues. However, this step is not necessary, the image will be pulled automatically when running the pipeline the first time.

```sh
docker pull nipengcsu/longmethyl:0.1
```

  - (2) Run longmethyl using `-profile docker`.

```sh
# activate nextflow environment
conda activate nextflow

# run longmethyl using cpu
nextflow run ~/tools/longmethyl -profile docker \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
# or, run longmethyl using GPU, set CUDA_VISIBLE_DEVICES and --gpu
CUDA_VISIBLE_DEVICES=0 nextflow run ~/tools/longmethyl -profile docker --gpu true \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
```

__Related issues__:

  1. For `No swap limit support`

```shell
# for Ubuntu

# (1) sudo, Edit the /etc/default/grub file. Add or edit the GRUB_CMDLINE_LINUX line 
# to add the following two key-value pairs
GRUB_CMDLINE_LINUX="cgroup_enable=memory swapaccount=1"

# (2) Update GRUB
sudo update-grub

# (3) Restart the machine
sudo reboot
```

Ref: [https://unix.stackexchange.com/questions/342735/docker-warning-no-swap-limit-support](https://unix.stackexchange.com/questions/342735/docker-warning-no-swap-limit-support)

  2. For `docker: Error response from daemon: could not select device driver "" with capabilities: [[gpu]].`

```shell
# for Ubuntu

distribution=$(. /etc/os-release;echo $ID$VERSION_ID) \
   && curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add - \
   && curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | sudo tee /etc/apt/sources.list.d/nvidia-docker.list

sudo apt-get update

sudo apt-get install -y nvidia-docker2

sudo systemctl restart docker
```

Ref: [https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)


  3. For `Failed to initialize NVML: Driver/library version mismatch`

Ref: [https://github.com/NVIDIA/nvidia-docker/issues/584](https://github.com/NVIDIA/nvidia-docker/issues/584)


### Option 3. Run with conda

  - (1) Install the conda environment named longmethyl (once for all).

```sh
# in a gpu machine, make sure there is already cuda10.0 driver in the machine
conda env create -f longmethyl/environment.yml
# or, in a cpu-only machine
conda env create -f longmethyl/environment-cpu.yml
```

  - (2) Install Guppy, since Guppy is not open-sourced, from [ONT community](https://nanoporetech.com/community) (once for all).


  - (3) Run longmethyl using `-profile conda` and the longmethyl environment.

```sh
# activate nextflow environment
conda activate nextflow

# run longmethyl
nextflow run ~/tools/longmethyl -profile conda \
    --conda_name /home/nipeng/tools/miniconda3/envs/longmethyl \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
# or, run longmethyl using GPU, set CUDA_VISIBLE_DEVICES
CUDA_VISIBLE_DEVICES=0 nextflow run ~/tools/longmethyl -profile conda \
    --conda_name /home/nipeng/tools/miniconda3/envs/longmethyl \
    --dsname test \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
```

## Acknowledgements
  - Some code were taken from [nanome](https://github.com/TheJacksonLaboratory/nanome) and [nf-core](https://github.com/nf-core).

developement: [nextflow_develop.md](docs/nextflow_develop.md)

## TODO
- add summmary
- test case with no basecall/resquiggle steps
- ~~dockerfile~~
- cpu settings (do not use task.cpus for all process)
- clean work dir
- ~~test with gpu~~ (singularity tested, still waiting to test docker+gpu, anyway)
- how to set a default deepsignal model

