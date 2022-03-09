#!/usr/bin/env nextflow
/*
========================================================================================
    longmethyl
========================================================================================
    Github : https://github.com/
----------------------------------------------------------------------------------------
*/

if( nextflow.version.matches(">= 20.07.1") ){
    nextflow.enable.dsl=2
} else {
    // Support lower version of nextflow
    nextflow.preview.dsl=2
}


/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

def helpMessage() {
    log.info"""
    longmethyl - Nextflow PIPELINE (v$workflow.manifest.version)
    =================================
    Usage:
    The typical command is as follows:

    nextflow run ~/tools/longmethyl -profile conda --conda_name /home/nipeng/tools/miniconda3/envs/longmethyl --genome GCF_000146045.2_R64_genomic.fna --input fast5s.al.demo/ --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz

    Mandatory arguments:
      --input       Input path for raw fast5 files (folders, tar/tar.gz files)
      --genome      Genome reference name ('hg38', 'ecoli', or 'hg38_chr22') or a directory, the directory must contain only one .fasta file with .fasta.fai index file. Default is hg38
      --dsname      Dataset/analysis name

    General options:
      --outdir      Output dir, default is 'results'
      --chrSet      Chromosomes used in analysis, default is chr1-22, X and Y, for human. For E. coli data, it is default as 'NC_000913.3'. For other reference genome, please specify each chromosome with space seperated.
      --cleanAnalyses   If clean old basecalling info in fast5 files
      --cleanup     If clean work dir after complete, default is false

    Tools specific options:
      --guppyDir        Guppy installation local directory, used only for conda environment
      --GUPPY_BASECALL_MODEL    Guppy basecalling model, default is 'dna_r9.4.1_450bps_hac.cfg'
      --deepsignalDir   DeepSignal model dir, default will get online
      --tomboResquiggleOptions  Tombo resquiggle options for super long/damaged sequencing, set to '--signal-length-range 0 500000  --sequence-length-range 0 50000'

    Running environment options:
      --docker_name     Docker name used for pipeline, default is 'liuyangzzu/nanome:latest'
      --singularity_name    Singularity name used for pipeline, default is 'docker://liuyangzzu/nanome:latest'
      --singularity_cache   Singularity cache dir, default is 'local_singularity_cache'
      --conda_name      Conda name used for pipeline, default is 'nanome'
      --conda_base_dir  Conda base directory, default is '/opt/conda'
      --conda_cache     Conda cache dir, default is 'local_conda_cache'

    -profile options:
      Use this parameter to choose a predefined configuration profile. Profiles can give configuration presets for different compute environments.

      test      A bundle of input params for ecoli test
      docker    A generic configuration profile to be used with Docker, pulls software from Docker Hub: liuyangzzu/nanome:latest
      singularity   A generic configuration profile to be used with Singularity, pulls software from: docker://liuyangzzu/nanome:latest
      conda     Please only use conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity. Check our GitHub for how to install local conda environment

    """.stripIndent()
}


// check input =====================================================================
// Show help message
if (params.help){
    helpMessage()
    exit 0
}

if (!params.genome){
    exit 1, "--genome option not specified!"
}

genome_map = params.genomes

if (params.genome && genome_map[params.genome]) { genome_path = genome_map[params.genome] }
else {  genome_path = params.genome }


// infer dataType, chrSet based on reference genome name, hg - human, ecoli - ecoli, otherwise is other reference genome
if (params.genome.contains('hg') || params.genome.contains('GRCh38') || params.genome.contains('GRCh37') || (params.dataType && params.dataType == 'human')) {
    dataType = "human"
    if (!params.chrSet) {
        // default for human, if false or 'false' (string), using '  '
        chrSet = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'
    } else {
        chrSet = params.chrSet
    }
} else if (params.genome.contains('ecoli') || (params.dataType && params.dataType == 'ecoli')) {
    dataType = "ecoli"
    if (!params.chrSet) {
        // default for ecoli
        chrSet = 'NC_000913.3'
    } else {
        chrSet = params.chrSet
    }
} else {
    // default will not found name, use other
    if (!params.dataType) { dataType = 'other' } else { dataType = params.dataType }
    if (!params.chrSet) {
        // No default value for other reference genome
        if (!file(params.genome+".contig_names.txt").exists()){
            // exit 1, "Missing --chrSet option for other reference genome, please specify chromosomes used in reference genome [${params.genome}], or use utils/extract_contig_names_from_fasta.py to create one"
        } else {
            Channel.fromPath( params.genome+".contig_names.txt", type: 'file', checkIfExists: true )
            .first(String)
            .set{chrSet}
        }
        
    }
    else { 
        chrSet = params.chrSet 
    }
}


// Collect all folders of fast5 files, and send into Channels for pipelines
if (params.input.endsWith(".filelist.txt")) {
    // list of files in filelist.txt
    Channel.fromPath( params.input, checkIfExists: true )
        .splitCsv(header: false)
        .map {
            if (!file(it[0]).exists())  {
                log.warn "File not exists: ${it[0]}, check file list: ${params.input}"
            } else {
                return file(it[0])
            }
        }
        .set{ fast5_tar_ch }
} else if (params.input.contains('*') || params.input.contains('?')) {
    // match all files in the folder, note: input must use '', prevent expand in advance
    // such as --input '/fastscratch/liuya/nanome/NA12878/NA12878_CHR22/input_chr22/*'
    Channel.fromPath(params.input, type: 'any', checkIfExists: true)
        .set{ fast5_tar_ch }
} else {
    // For single file/wildcard matched files
    Channel.fromPath( params.input, checkIfExists: true ).set{ fast5_tar_ch }
}


// set utils/src dirs 
projectDir = workflow.projectDir
ch_utils = Channel.fromPath("${projectDir}/utils",  type: 'dir', followLinks: false)
ch_src   = Channel.fromPath("${projectDir}/src",  type: 'dir', followLinks: false)


// TODO: set summary
def summary = [:]
summary['input']            = params.input


// Reference genome
def referenceGenome = 'reference_genome/ref.fasta'

// Check all tools work well; from nanome
process EnvCheck {
    tag "envcheck"
    errorStrategy 'terminate'

    publishDir "${params.outdir}/${params.dsname}-ds",
        mode: "copy", pattern: "tools_version_table.tsv", overwrite: true

    input:
    path utils
    path deepsignalDir
    path reference_genome

    output:
    path "reference_genome",                emit: reference_genome, optional: true
    path "${params.DEEPSIGNAL_MODEL_DIR}",  emit: deepsignal_model, optional: true

    script:
    """
    date; hostname; pwd
    echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"

    ## Validate nanome container/environment is correct
    # bash utils/validate_nanome_container.sh  tools_version_table.tsv

    ## Untar and prepare deepsignal model
    if [ ${params.runDeepSignal} == true ]; then
        if [ ${deepsignalDir} == *.tar.gz ] ; then
            ## Get DeepSignal Model online
            tar -xzf ${deepsignalDir}
        elif [[ ${deepsignalDir} != ${params.DEEPSIGNAL_MODEL_DIR} && -d ${deepsignalDir} ]] ; then
            ## rename it to deepsignal default dir name
            cp -a ${deepsignalDir}  ${params.DEEPSIGNAL_MODEL_DIR}
        fi
        ## Check DeepSignal model
        ls -lh ${params.DEEPSIGNAL_MODEL_DIR}/
    fi

    if [[ ${params.runBasecall} == true || ${params.runDeepSignal} == true ]]; then
        ## Get dir for reference_genome
        mkdir -p reference_genome
        find_dir="\$PWD/reference_genome"
        if [[ ${reference_genome} == *.tar.gz && -f ${reference_genome}  ]] ; then
            tar -xzf ${reference_genome} -C reference_genome
        elif [[ ${reference_genome} == *.tar && -f ${reference_genome} ]] ; then
            tar -xf ${reference_genome} -C reference_genome
        elif [[ -d ${reference_genome} ]] ; then
            ## for folder, use ln, note this is a symbolic link to a folder
            find_dir=\$( readlink -f ${reference_genome} )
        elif [[ -e ${reference_genome} ]] ; then
            cp ${reference_genome} \${find_dir}/${reference_genome}.ori.fasta
        else
            echo "### ERROR: not recognized reference_genome=${reference_genome}"
            exit -1
        fi

        find \${find_dir} -name '*.fasta*' | \
             parallel -j0 -v  'fn={/} ; ln -s -f  {}   reference_genome/\${fn/*.fasta/ref.fasta}'
        ## find \${find_dir} -name '*.sizes' | \
        ##         parallel -j1 -v ln -s -f {} reference_genome/chrom.sizes

        ls -lh reference_genome/
    fi

    echo "### Check env"
    echo "cpus=$task.cpus"
    echo "referenceGenome=${referenceGenome}"
    echo "### Check env DONE"
    """
}


// Untar of subfolders named 'M1', ..., 'M10', etc.
process Untar {
    tag "${fast5_tar.baseName}"

    label 'process_medium'

    input:
    path fast5_tar

    output:
    path "${fast5_tar.baseName}.untar", emit:untar,  optional: true

    script:
    cores = task.cpus
    if (!params.runBasecall) { // perform basecall
        """
        date; hostname; pwd
        echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"

        ## Extract input files tar/tar.gz/folder
        mkdir -p untarTempDir
        if [[ ${fast5_tar} == *.tar && -f ${fast5_tar} ]] ; then
            ### deal with tar
            tar -xf ${fast5_tar} -C untarTempDir
        elif [[ ${fast5_tar} == *.tar.gz && -f ${fast5_tar} ]] ; then
            ### deal with tar.gz
            tar -xzf ${fast5_tar} -C untarTempDir
        elif [[ -d ${fast5_tar} ]]; then
            ## Copy files, do not change original files such as old analyses data
            find ${fast5_tar}/ -name '*.fast5' | \
                parallel -j$cores  cp {} untarTempDir/
        else
            echo "### Untar error for input=${fast5_tar}"
        fi

        ## Move fast5 raw/basecalled files into XXX.untar folder
        mkdir -p ${fast5_tar.baseName}.untar

        find untarTempDir -name "*.fast5" -type f | \
            parallel -j$cores  mv {}  ${fast5_tar.baseName}.untar/

        ## Clean temp files
        rm -rf untarTempDir

        ## Clean old basecalled analyses in input fast5 files
        if [[ "${params.cleanAnalyses}" == true ]] ; then
            echo "### Start cleaning old analysis"
            ## python -c 'import h5py; print(h5py.version.info)'
            clean_old_basecall_in_fast5.py \
                -i ${fast5_tar.baseName}.untar --is-indir --verbose\
                --processor $cores
        fi

        totalFiles=\$( find ${fast5_tar.baseName}.untar -name "*.fast5" -type f | wc -l )
        echo "### Total fast5 input files:\${totalFiles}"
        if (( totalFiles==0 )); then
            echo "### no fast5 files at ${fast5_tar.baseName}.untar, skip this job"
            rm -rf ${fast5_tar.baseName}.untar
        fi
        echo "### Untar DONE"
        """
    } else {
        """
        date; hostname; pwd

        ## Extract input files tar/tar.gz/folder
        mkdir -p untarTempDir
        if [[ ${fast5_tar} == *.tar && -f ${fast5_tar} ]] ; then
            ### deal with tar
            tar -xf ${fast5_tar} -C untarTempDir
        elif [[ ${fast5_tar} == *.tar.gz && -f ${fast5_tar} ]] ; then
            ### deal with tar.gz
            tar -xzf ${fast5_tar} -C untarTempDir
        elif [[ -d ${fast5_tar} ]] ; then
            ## user provide basecalled input dir, just cp them
            mkdir -p untarTempDir/test
            cp -rf ${fast5_tar}/*   untarTempDir/test/
        else
            echo "### Untar error for input=${fast5_tar}"
        fi

        ## Move fast5 raw/basecalled files into XXX.untar folder
        mkdir -p ${fast5_tar.baseName}.untar
        ## Keep the directory structure for basecalled input
        mv untarTempDir/*/*   ${fast5_tar.baseName}.untar/

        ## Clean temp files
        rm -rf untarTempDir

        totalFiles=\$( find ${fast5_tar.baseName}.untar -name "*.fast5" -type f | wc -l )
        echo "### Total fast5 input files:\${totalFiles}"
        if (( totalFiles==0 )); then
            echo "### no fast5 files at ${fast5_tar.baseName}.untar, skip this job"
            rm -rf ${fast5_tar.baseName}.untar
        fi
        echo "### Untar DONE"
        """
    }
}


// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
    tag "${fast5_dir.baseName}"

    label 'process_medium'

    input:
    path fast5_dir

    output:
    path "${fast5_dir.baseName}.basecalled",    emit: basecall

    when:
    params.runBasecall

    script:
    cores = task.cpus

    """
    date; hostname; pwd

    echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"
    if [[ "\${CUDA_VISIBLE_DEVICES:-}" == "" ]] ; then
        echo "Detect no GPU, using CPU commandType"
        commandType='cpu'
        gpuOptions=" "
    else
        echo "Detect GPU, using GPU commandType"
        commandType='gpu'
        gpuOptions="-x auto"
    fi

    which guppy_basecaller
    guppy_basecaller -v
    mkdir -p ${fast5_dir.baseName}.basecalled

    if [[ ${params.runBasecall} == true ]] ; then
        ## CPU/GPU version command
        guppy_basecaller --input_path ${fast5_dir} \
            --save_path "${fast5_dir.baseName}.basecalled" \
            --config ${params.GUPPY_BASECALL_MODEL} \
            --num_callers ${cores} \
            --fast5_out --compress_fastq\
            --verbose_logs  \${gpuOptions} &>> ${params.dsname}.${fast5_dir.baseName}.Basecall.run.log
    else
        ## Just use user's basecalled input
        cp -rf ${fast5_dir}/*   ${fast5_dir.baseName}.basecalled/
    fi

    ## Combine fastq
    touch "${fast5_dir.baseName}.basecalled"/batch_basecall_combine_fq_${fast5_dir.baseName}.fq.gz

    ## Below is compatable with both Guppy v4.2.2 (old) and newest directory structures
    if [[ -d ${fast5_dir.baseName}.basecalled/pass && -d ${fast5_dir.baseName}.basecalled/fail ]]; then
        find "${fast5_dir.baseName}.basecalled/" "${fast5_dir.baseName}.basecalled/pass/"\
            "${fast5_dir.baseName}.basecalled/fail/" -maxdepth 1 -name '*.fastq.gz' -type f\
            -print0 2>/dev/null | \
            while read -d \$'\0' file ; do
                cat \$file >> \
                    "${fast5_dir.baseName}.basecalled"/batch_basecall_combine_fq_${fast5_dir.baseName}.fq.gz
            done
    else
        find "${fast5_dir.baseName}.basecalled/" -maxdepth 1 -name '*.fastq.gz' -type f\
            -print0 2>/dev/null | \
            while read -d \$'\0' file ; do
                cat \$file >> \
                    "${fast5_dir.baseName}.basecalled"/batch_basecall_combine_fq_${fast5_dir.baseName}.fq.gz
            done
    fi
    echo "### Combine fastq.gz DONE"

    ## Remove fastq.gz
    if [[ -d ${fast5_dir.baseName}.basecalled/pass && -d ${fast5_dir.baseName}.basecalled/fail ]]; then
        find "${fast5_dir.baseName}.basecalled/"   "${fast5_dir.baseName}.basecalled/pass/"\
            "${fast5_dir.baseName}.basecalled/fail/" -maxdepth 1 -name '*.fastq.gz' -type f 2>/dev/null |\
            parallel -j${cores} 'rm -f {}'
    else
        find "${fast5_dir.baseName}.basecalled/" -maxdepth 1 -name '*.fastq.gz' -type f 2>/dev/null |\
            parallel -j${cores} 'rm -f {}'
    fi

    ## After basecall, rename and publish summary filenames, summary may also be used by resquiggle
    mv ${fast5_dir.baseName}.basecalled/sequencing_summary.txt \
        ${fast5_dir.baseName}.basecalled/${fast5_dir.baseName}-sequencing_summary.txt

    echo "### Basecalled by Guppy DONE"
    """
}


// Resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
    tag "${basecallIndir.baseName}"

    label 'process_high'

    input:
    path    basecallIndir
    each    path(reference_genome)

    output:
    path "${basecallIndir.baseName}.resquiggle",    emit: resquiggle

    when:
    params.runDeepSignal

    script:
    cores = task.cpus
    samtools_cores = task.cpus
    resquiggle_cores = (task.cpus).intValue()
    
    """
    rm -rf ${basecallIndir.baseName}.resquiggle
    mkdir -p ${basecallIndir.baseName}.resquiggle/workspace
    cp -f ${basecallIndir}/batch_basecall_combine_fq_*.fq.gz  \
        ${basecallIndir.baseName}.resquiggle/
    cp -f ${basecallIndir}/${basecallIndir.baseName}-sequencing_summary.txt  \
        ${basecallIndir.baseName}.resquiggle/
    find ${basecallIndir}/workspace -name '*.fast5' -type f| \
        parallel -j${cores}  \
        'cp {}   ${basecallIndir.baseName}.resquiggle/workspace/'
    echo '### Duplicate from basecall DONE'
    gunzip ${basecallIndir.baseName}.resquiggle/batch_basecall_combine_fq_*.fq.gz
    tombo preprocess annotate_raw_with_fastqs\
        --fast5-basedir ${basecallIndir.baseName}.resquiggle/workspace\
        --fastq-filenames ${basecallIndir.baseName}.resquiggle/batch_basecall_combine_fq_*.fq\
        --sequencing-summary-filenames ${basecallIndir.baseName}.resquiggle/${basecallIndir.baseName}-sequencing_summary.txt\
        --basecall-group ${params.BasecallGroupName}\
        --basecall-subgroup ${params.BasecallSubGroupName}\
        --overwrite --processes  ${samtools_cores} \
        &>> ${params.dsname}.${basecallIndir.baseName}.Resquiggle.run.log
    echo '### tombo preprocess DONE'
    tombo resquiggle\
        --processes ${resquiggle_cores} \
        --corrected-group ${params.ResquiggleCorrectedGroup} \
        --basecall-group ${params.BasecallGroupName} \
        --basecall-subgroup ${params.BasecallSubGroupName}\
        --ignore-read-locks ${params.tomboResquiggleOptions ? params.tomboResquiggleOptions : ''} \
        --overwrite \
        ${basecallIndir.baseName}.resquiggle/workspace \
        ${referenceGenome} &>> ${params.dsname}.${basecallIndir.baseName}.Resquiggle.run.log  
    echo '### tombo resquiggle DONE'
    """
}


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
    tag "${indir.baseName}"

    label 'process_high'

    publishDir "${params.outdir}/${params.dsname}_intermediate/deepsignal",
        mode: "copy",
        enabled: params.outputIntermediate

    input:
    path indir
    each path(reference_genome)
    each path(deepsignal_model_dir)

    output:
    path "${params.dsname}_deepsignal_batch_${indir.baseName}.*.gz",    emit: deepsignal_out

    when:
    params.runDeepSignal

    script:
    cores = task.cpus
    """
    DeepSignalModelBaseDir="."
    ### commandType='cpu'
    outFile="${params.dsname}_deepsignal_batch_${indir.baseName}.tsv"

    if [[ "\${CUDA_VISIBLE_DEVICES:-}" == "" ]] ; then
        echo "Detect no GPU, using CPU commandType"
        commandType='cpu'
    else
        echo "Detect GPU, using GPU commandType"
        commandType='gpu'
    fi

    if [[ \${commandType} == "cpu" ]]; then
        ## CPU version command
        deepsignal call_mods \
            --input_path ${indir}/workspace \
            --model_path "\${DeepSignalModelBaseDir}/${params.DEEPSIGNAL_MODEL_DIR}/${params.DEEPSIGNAL_MODEL}" \
            --result_file \${outFile} \
            --reference_path ${referenceGenome} \
            --corrected_group ${params.ResquiggleCorrectedGroup} \
            --nproc $cores \
            --is_gpu no   &>> ${params.dsname}.${indir.baseName}.DeepSignal.run.log
    elif [[ \${commandType} == "gpu" ]]; then
        ## GPU version command
        deepsignal call_mods \
            --input_path ${indir}/workspace \
            --model_path "\${DeepSignalModelBaseDir}/${params.DEEPSIGNAL_MODEL_DIR}/${params.DEEPSIGNAL_MODEL}" \
            --result_file \${outFile} \
            --reference_path ${referenceGenome} \
            --corrected_group ${params.ResquiggleCorrectedGroup} \
            --nproc $cores \
            --is_gpu yes  &>> ${params.dsname}.${indir.baseName}.DeepSignal.run.log
    else
        echo "### error value for commandType=\${commandType}"
        exit 255
    fi

    gzip -f \${outFile}
    echo "### DeepSignal methylation DONE"
    """
}


// Combine DeepSignal runs' all results together
process DeepSignalFreq {
    tag "${params.dsname}"

    publishDir "${params.outdir}/${params.dsname}-ds",
        mode: "copy",
        pattern: "${params.dsname}_deepsignal_per_read_combine.*.gz",
        enabled: params.outputRaw

    publishDir "${params.outdir}/${params.dsname}-ds",
        mode: "copy", pattern: "${params.dsname}_deepsignal_sitemods_freq.bed.gz"

    input:
    path x
    path ch_utils
    path ch_src

    output:
    path "${params.dsname}_deepsignal_per_read_combine.*.gz",   emit: deepsignal_combine_out
    path "${params.dsname}_deepsignal_sitemods_freq.bed.gz",   emit: site_unify

    when:
    x.size() >= 1 && params.runCallfreq

    """
    touch ${params.dsname}_deepsignal_per_read_combine.tsv.gz
    cat ${x} > ${params.dsname}_deepsignal_per_read_combine.tsv.gz

    if [[ ${params.deduplicate} == true ]] ; then
        echo "### Deduplicate for read-level outputs"
        ## sort order: Chr, Start, (End), ID, Strand
        zcat ${params.dsname}_deepsignal_per_read_combine.tsv.gz |\
            sort -V -u -k1,1 -k2,2n -k5,5 -k3,3 |\
            gzip -f > ${params.dsname}_deepsignal_per_read_combine.sort.tsv.gz
        rm ${params.dsname}_deepsignal_per_read_combine.tsv.gz &&\
            mv ${params.dsname}_deepsignal_per_read_combine.sort.tsv.gz  \
                ${params.dsname}_deepsignal_per_read_combine.tsv.gz
    fi

    ## Unify format output
    bash utils/call_freq_deepsignal.sh \
        ${params.dsname} deepsignal \
        ${params.dsname}_deepsignal_per_read_combine.tsv.gz \
        ${params.sort  ? true : false}
    echo "### DeepSignal combine DONE"
    """
}


/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    if ( !file(params.genome).exists() )
        exit 1, "genome reference path does not exist, check params: --genome ${params.genome}"

    genome_ch = Channel.fromPath(genome_path, type: 'any', checkIfExists: true)

    if (! params.runDeepSignal) {
        // use null placeholder
        deepsignalDir = Channel.fromPath("${projectDir}/utils/null2", type: 'any', checkIfExists: true)
    }
    else {
        // User provide the dir
        if ( !file(params.deepsignalDir.toString()).exists() )
            exit 1, "deepsignalDir does not exist, check params: --deepsignalDir ${params.deepsignalDir}"
        deepsignalDir = Channel.fromPath(params.deepsignalDir, type: 'any', checkIfExists: true)
    }

    EnvCheck(ch_utils, deepsignalDir, genome_ch)
    Untar(fast5_tar_ch)
    Basecall(Untar.out.untar)

    // Resquiggle running if use Tombo or DeepSignal
    Resquiggle(Basecall.out.basecall, EnvCheck.out.reference_genome)

    if (params.runDeepSignal) {
        DeepSignal(Resquiggle.out.resquiggle, EnvCheck.out.reference_genome, EnvCheck.out.deepsignal_model)
        comb_deepsignal = DeepSignalFreq(DeepSignal.out.deepsignal_out.collect(), ch_utils, ch_src)
        s3 = comb_deepsignal.site_unify
    } else {
        s3 = Channel.empty()
    }
    
}

/*
========================================================================================
    THE END
========================================================================================
*/
