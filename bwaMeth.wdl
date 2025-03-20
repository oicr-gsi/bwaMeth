version 1.0

struct GenomeResources {
    String indexModule 
    String index
    String genomeModule
    String fasta
}

struct fastqInputs {
    File fastqR1
    File? fastqR2
    String readGroup
    Int numChunk
    Int? numReads
}

workflow bwaMeth {
    input {
        Array[fastqInputs] inputGroups
        String outputFileNamePrefix     
        String reference
    }

    parameter_meta {
        inputGroups: "Array of fastq inputs and parameters"
        outputFileNamePrefix: "Prefix for output files"
        reference: "The genome reference build. For example: hg19, hg38"
    }

    Map[String, GenomeResources] resources = {
        "hg38": {
            "indexModule": "hg38-bwa-meth-index/p12-2022-10-17",
            "index": "$HG38_BWA_METH_INDEX_ROOT/hg38_random.fa",
            "genomeModule": "hg38-em-seq/p12-2022-10-17",
            "fasta": "$HG38_EM_SEQ_ROOT/hg38_random.fa"
        }
    }

    GenomeResources ref = resources[reference]

    scatter (ig in inputGroups){
        if (ig.numChunk > 1) {
            call countChunkSize {
                input:
                fastqR1 = ig.fastqR1,
                numChunk = ig.numChunk,
                numReads = ig.numReads
            }
        
            call slicer as slicerR1 { 
                input: 
                fastqR = ig.fastqR1,
                chunkSize = countChunkSize.chunkSize
            }
            if (defined(ig.fastqR2)) {
                File fastqR2_ = select_all([ig.fastqR2])[0]
                call slicer as slicerR2 {
                    input:
                    fastqR = fastqR2_,
                    chunkSize = countChunkSize.chunkSize
                }
            }
        }

        Array[File] fastq1 = select_first([slicerR1.chunkFastq, [ig.fastqR1]])

        if(defined(ig.fastqR2)) {
        Array[File?] fastq2 = select_first([slicerR2.chunkFastq, [ig.fastqR2]])
        Array[Pair[File,File?]] pairedFastqs = zip(fastq1,fastq2)
        }

        if(!defined(ig.fastqR2)) {
        Array[Pair[File,File?]] singleFastqs = cross(fastq1,[ig.fastqR2])
        }

        Array[Pair[File,File?]] fastqPairs = select_first([pairedFastqs, singleFastqs])

        scatter (p in fastqPairs) {
            call trimAndAlign  { 
                    input: 
                    read1 =  p.left,
                    read2 = if (defined(ig.fastqR2)) then p.right else ig.fastqR2,
                    bwaReadGroup = ig.readGroup,
                    bwaIndex = ref.index,
                    modules = "fastp/0.23.2 bwa-meth/0.2.5 ~{ref.indexModule}"
            }    
        }
        call mergeBams {
            input:
            bams = trimAndAlign.bam
        }
    }
    
    call mergeAandMarkDuplicates {
        input:
        bams = mergeBams.mergedBam,
        outputFileNamePrefix = outputFileNamePrefix,
        reference_genome = ref.index,
        modules = "sambamba/0.8.2 samtools/1.15  picard/2.21.2 methylseq-mark-nonconverted-reads/1.2 ~{ref.indexModule}"
    }

    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "Workflow to run bwa-meth, the fast aligner for EM-seq/BS-Seq reads. Prior to alignment, adatper trimming and quality filtering are performed. Readgroup information to be injected into the bam header needs to be provided.  The workflow can also split the input data into a requested number of chunks, align each separately then merge the separate alignments into a single bam file.  This decreases the workflow run time. Final bam file also applied markDuplicates."
        dependencies: [
        {
            name: "fastp/0.23.2",
            url: "https://github.com/OpenGene/fastp"
        },
        {
            name: "bwa-meth/0.2.5",
            url: "https://github.com/brentp/bwa-meth"
        },
        {
            name: "slicer/0.3.0",
            url: "https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz"
        },
        { 
          name: "samtools/1.15",
          url: "https://github.com/samtools/samtools/releases/"
        },
        { 
          name: "gsi hg38 modules : hg38-em-seq/p12-2022-10-17",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        },
        {
          name: "gsi modules : hg38-bwa-meth-index/p12-2022-10-17",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        },
        { name: "picard/2.21.2",
          url: "https://broadinstitute.github.io/picard/"
        },
        {
            name: "python/3.7",
            url: "https://www.python.org"
        },
        {   name: "methylseq-mark-nonconverted-reads/1.2",
            url: "https://github.com/nebiolabs/mark-nonconverted-reads"
        }
      ]
      output_meta: {
        bwaMethBam: {
            description: "Output Alignment BAM file, merged, dedplicated and marked nonconverted reads",
            vidarr_label: "bwaMethBam"
        },
        bwaMethBamIndex: {
            description: "Index of the Output Alignment file, merged, deduplicated file, and marked nonconverted reads",
            vidarr_label: "bwaMethIndex"
        },
        nonConvertedReads: {
            description: "Statistics of nonconverted reads",
            vidarr_label: "nonConvertedReads"
        }
      }
    }

    output {
        File bwaMethBam = mergeAandMarkDuplicates.outputMergedBam
        File bwaMethBamIndex = mergeAandMarkDuplicates.outputMergedBai
        File nonConvertedReads = mergeAandMarkDuplicates.nonconverted_reads
    }
}


task countChunkSize{
    input {
        File fastqR1
        Int numChunk
        Int? numReads
        String modules = "python/3.7"
        Int jobMemory = 16
        Int timeout = 48
    }
    
    parameter_meta {
        fastqR1: "Fastq file for read 1"
        numChunk: "Number of chunks to split fastq file"
        numReads: "Number of reads"
        modules: "Required environment modules"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
    
    command <<<
        set -euo pipefail

        if [ -z "~{numReads}" ]; then
            totalLines=$(zcat ~{fastqR1} | wc -l)
        else totalLines=$((~{numReads}*4))
        fi
        
        python3 -c "from math import ceil; print (int(ceil(($totalLines/4.0)/~{numChunk})*4))"
    >>>
    
    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }
    
    output {
        String chunkSize =read_string(stdout())
    }

    meta {
    output_meta: {
      chunkSize: "output number of lines per chunk"
    }
    }    
   
}

task slicer {
    input {
        File fastqR         
        String chunkSize
        String modules = "slicer/0.3.0"
        Int jobMemory = 16
        Int timeout = 48
    }
    
    parameter_meta {
        fastqR: "Fastq file"
        chunkSize: "Number of lines per chunk"
        modules: "Required environment modules"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
    
    command <<<
        set -euo pipefail
        slicer -i ~{fastqR} -l ~{chunkSize} --gzip 
    >>>
    
    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    } 
    
    output {
        Array[File] chunkFastq = glob("*.fastq.gz")
    }

    meta {
        output_meta: {
            chunkFastq: "output fastq chunks"
        }
    } 
  
}

task trimAndAlign {
    input {
        File read1
        File? read2

        String bwaReadGroup
        String bwaIndex

        Boolean fastpDisableQualityFiltering = false
        Int? fastpQualifiedQualityPhred
        Int? fastpUnqualifiedPercentLimit
        Int? fastpNBaseLimit

        Boolean fastpDisableLengthFiltering = false
        Int? fastpLengthRequired

        Boolean fastpDisableAdapterTrimming = false

        Boolean fastpDisableTrimPolyG = false

        Int timeout = 48
        Int memory = 32
        Int threads = 8
        String modules
    }

    parameter_meta {
        read1: "Read 1 FastQ file"
        read2: "Read 2 FastQ file"
        bwaReadGroup: "Read group that will populate the `@RG` BAM flag"
        bwaIndex: "The FastA in the directory that contains the bwa index files"
        fastpDisableQualityFiltering: "Disable fastp quality filtering"
        fastpQualifiedQualityPhred: "The quality value that a base is considered qualified (default >=Q15)"
        fastpUnqualifiedPercentLimit: "How many percents of bases are allowed to be unqualified (default 40%)"
        fastpNBaseLimit: "How many N can a read have before being discarded (default 5)"
        fastpDisableLengthFiltering: "Disable filtering reads below a certain length"
        fastpLengthRequired: "Reads shorter than length_required will be discarded (default 15)"
        fastpDisableAdapterTrimming: "Disable all adapter trimming"
        fastpDisableTrimPolyG: "Disable triming polyG at the end of the read"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    String fastpQ = if fastpDisableQualityFiltering then "-Q" else ""
    String fastpq = if defined(fastpQualifiedQualityPhred) then "-q ~{fastpQualifiedQualityPhred}" else ""
    String fastpu = if defined(fastpUnqualifiedPercentLimit) then "-u ~{fastpUnqualifiedPercentLimit}" else ""
    String fastpn = if defined(fastpNBaseLimit) then "-n ~{fastpNBaseLimit}" else ""

    String fastpL = if fastpDisableLengthFiltering then "-L" else ""
    String fastpl = if defined(fastpNBaseLimit) then "-l ~{fastpNBaseLimit}" else ""

    String fastpA = if fastpDisableAdapterTrimming then "-A" else ""

    String fastpG = if fastpDisableTrimPolyG then "-G" else ""

    command <<<
        set -euo pipefail
        fastp \
            --stdout --thread ~{threads} \
            ~{fastpQ} ~{fastpq} ~{fastpu} ~{fastpn} ~{fastpL} ~{fastpl} ~{fastpA} ~{fastpG} \
            -i ~{read1} -I ~{read2} \
        | bwameth.py -p --threads ~{threads} --read-group ~{bwaReadGroup} --reference ~{bwaIndex} /dev/stdin \
        | samtools sort -o output.bam -@ ~{threads} -
    >>>

    output {
        File fastpReport = "fastp.json"
        File bam = "output.bam"
    }

    meta {
        output_meta: {
            fastpReport: "The json report file produced by fastp",
            bam: "The bam file produced by the trimmed FastQ files fed to bwa-meth"
        }
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task mergeBams{
    input {
        Array[File] bams
        Int jobMemory = 32
        String modules = "picard/2.21.2"
        Int timeout = 12
    }
    parameter_meta {
        bams:  "Input bam files"
        jobMemory: "Memory allocated indexing job"
        modules:   "Required environment modules"
        timeout:   "Hours before task timeout"    
    }

    command <<<
        set -euo pipefail

        export JAVA_OPTS="-Xmx$(echo "scale=0; ~{jobMemory} * 0.8 / 1" | bc)G"
        java -jar ${PICARD_ROOT}/picard.jar \
        MergeSamFiles \
        I=~{sep=" I=" bams} \
        O=mergedChunks.bam \
        USE_THREADING=true \
        SORT_ORDER=coordinate
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File mergedBam = "mergedChunks.bam"
    }
}

task mergeAandMarkDuplicates{
    input {
        Array[File] bams
        String outputFileNamePrefix
        String reference_genome
        Int opticalDistance = 100
        Int jobMemory = 64
        String modules
        Int timeout = 72
    }
    parameter_meta {
        bams:  "Input bam files"
        outputFileNamePrefix: "Prefix for output file"
        reference_genome: "the reference genome fasta"
        opticalDistance: "For MarkDuplicates. The maximum offset between two duplicate clusters in order to consider them optical duplicates. 100 is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more appropriate."
        jobMemory: "Memory allocated indexing job"
        modules:   "Required environment modules"
        timeout:   "Hours before task timeout"    
    }
    String tmpDir = "tmp/"

    command <<<
        set -euo pipefail
        mkdir -p ~{tmpDir}

        export JAVA_OPTS="-Xmx$(echo "scale=0; ~{jobMemory} * 0.8 / 1" | bc)G"
        java -jar ${PICARD_ROOT}/picard.jar \
        MergeSamFiles \
        I=~{sep=" I=" bams} \
        O=~{outputFileNamePrefix}.merged.bam \
        USE_THREADING=true \
        SORT_ORDER=coordinate

        java -jar ${PICARD_ROOT}/picard.jar \
        MarkDuplicates \
        I=~{outputFileNamePrefix}.merged.bam \
        O=~{outputFileNamePrefix}.merged.deduped.bam \
        METRICS_FILE=~{outputFileNamePrefix}.markDuplicates.txt \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{opticalDistance} \
        CREATE_INDEX=true \
        ASSUME_SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=SILENT

        python3 $METHYLSEQ_MARK_NONCONVERTED_READS_ROOT/bin/mark-nonconverted-reads.py --reference ~{reference_genome} --bam ~{outputFileNamePrefix}.merged.deduped.bam 2> "~{outputFileNamePrefix}.merged.deduped.nonconverted.tsv" \
        | samtools view -u /dev/stdin \
        | sambamba sort  --tmpdir=~{tmpDir}  -o "~{outputFileNamePrefix}.marknonconverted.merged.deduped.bam" /dev/stdin
        sambamba index ~{outputFileNamePrefix}.marknonconverted.merged.deduped.bam
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File outputMergedBam = "~{outputFileNamePrefix}.marknonconverted.merged.deduped.bam"
        File outputMergedBai = "~{outputFileNamePrefix}.marknonconverted.merged.deduped.bam.bai"
        File nonconverted_reads = "~{outputFileNamePrefix}.merged.deduped.nonconverted.tsv"
    }

    meta {
        output_meta: {
            outputMergedBam: "output merged bam aligned to genome"
        }
    }       
}

