# longmethyl
A demo nextflow pipeline of methylation detection using long reads

<p>&nbsp;&nbsp;</p>

## Installation

  - (1) Install conda from [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) if neeeded.


  - (2) Install nextflow.

Create an environment containing nextflow/install nextflow:
```sh
# create a new environment
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

<p>&nbsp;&nbsp;</p>

## Usage

### Option 1. Run with docker


<p>&nbsp;</p>

### Option 2. Run with singularity



<p>&nbsp;&nbsp;</p>

### Option 3. Run with conda

  - (1) Install conda environment (once for all).

```sh
conda env create -f longmethyl/environment.yml
```

  - (2) Install Guppy, since Guppy is not open-sourced, from [ONT community](https://nanoporetech.com/community) (once for all).


  - (3) Run longmethyl in an environment containing nextflow.

```sh
## demo
# activate nextflow environment
conda activate nextflow
# run longmethyl
nextflow run ~/tools/longmethyl -profile conda --conda_name /home/nipeng/tools/miniconda3/envs/longmethyl \
    --genome GCF_000146045.2_R64_genomic.fna \
    --input fast5s.al.demo/ \
    --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
```

<p>&nbsp;</p>


## Acknowledgements
  - Some code were taken from [nanome](https://github.com/TheJacksonLaboratory/nanome) or [nf-core](https://github.com/nf-core).



For nextflow developers: [nextflow_develop.md](docs/nextflow_develop.md)


## TODO
- add summmary
- test case with no basecall/resquiggle steps
- dockerfile
- cpu settings
- clean work dir

