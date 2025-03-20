# bwaMeth

Workflow to run bwa-meth, the fast aligner for EM-seq/BS-Seq reads. Prior to alignment, adatper trimming and quality filtering are performed. Readgroup information to be injected into the bam header needs to be provided.  The workflow can also split the input data into a requested number of chunks, align each separately then merge the separate alignments into a single bam file.  This decreases the workflow run time. Final bam file also applied markDuplicates.

## Overview

## Dependencies

* [fastp 0.23.2](https://github.com/OpenGene/fastp)
* [bwa-meth 0.2.5](https://github.com/brentp/bwa-meth)
* [slicer 0.3.0](https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz)
* [samtools 1.15](https://github.com/samtools/samtools/releases/)
* [gsi hg38 modules : hg38-em-seq p12-2022-10-17](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi modules : hg38-bwa-meth-index p12-2022-10-17](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [picard 2.21.2](https://broadinstitute.github.io/picard/)
* [python 3.7](https://www.python.org)
* [methylseq-mark-nonconverted-reads 1.2](https://github.com/nebiolabs/mark-nonconverted-reads)


## Usage

### Cromwell
```
java -jar cromwell.jar run bwaMeth.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputGroups`|Array[fastqInputs]|Array of fastq inputs and parameters
`outputFileNamePrefix`|String|Prefix for output files
`reference`|String|The genome reference build. For example: hg19, hg38


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`countChunkSize.modules`|String|"python/3.7"|Required environment modules
`countChunkSize.jobMemory`|Int|16|Memory allocated for this job
`countChunkSize.timeout`|Int|48|Hours before task timeout
`slicerR1.modules`|String|"slicer/0.3.0"|Required environment modules
`slicerR1.jobMemory`|Int|16|Memory allocated for this job
`slicerR1.timeout`|Int|48|Hours before task timeout
`slicerR2.modules`|String|"slicer/0.3.0"|Required environment modules
`slicerR2.jobMemory`|Int|16|Memory allocated for this job
`slicerR2.timeout`|Int|48|Hours before task timeout
`trimAndAlign.fastpDisableQualityFiltering`|Boolean|false|Disable fastp quality filtering
`trimAndAlign.fastpQualifiedQualityPhred`|Int?|None|The quality value that a base is considered qualified (default >=Q15)
`trimAndAlign.fastpUnqualifiedPercentLimit`|Int?|None|How many percents of bases are allowed to be unqualified (default 40%)
`trimAndAlign.fastpNBaseLimit`|Int?|None|How many N can a read have before being discarded (default 5)
`trimAndAlign.fastpDisableLengthFiltering`|Boolean|false|Disable filtering reads below a certain length
`trimAndAlign.fastpLengthRequired`|Int?|None|Reads shorter than length_required will be discarded (default 15)
`trimAndAlign.fastpDisableAdapterTrimming`|Boolean|false|Disable all adapter trimming
`trimAndAlign.fastpDisableTrimPolyG`|Boolean|false|Disable triming polyG at the end of the read
`trimAndAlign.timeout`|Int|48|The hours until the task is killed
`trimAndAlign.memory`|Int|32|The GB of memory provided to the task
`trimAndAlign.threads`|Int|8|The number of threads the task has access to
`mergeBams.jobMemory`|Int|32|Memory allocated indexing job
`mergeBams.modules`|String|"picard/2.21.2"|Required environment modules
`mergeBams.timeout`|Int|12|Hours before task timeout
`mergeAandMarkDuplicates.opticalDistance`|Int|100|For MarkDuplicates. The maximum offset between two duplicate clusters in order to consider them optical duplicates. 100 is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more appropriate.
`mergeAandMarkDuplicates.jobMemory`|Int|64|Memory allocated indexing job
`mergeAandMarkDuplicates.timeout`|Int|72|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`bwaMethBam`|File|Output Alignment BAM file, merged, dedplicated and marked nonconverted reads|vidarr_label: bwaMethBam
`bwaMethBamIndex`|File|Index of the Output Alignment file, merged, deduplicated file, and marked nonconverted reads|vidarr_label: bwaMethIndex
`nonConvertedReads`|File|Statistics of nonconverted reads|vidarr_label: nonConvertedReads


## Commands
 This section lists command(s) run by bwaMeth workflow
 
 * Running bwaMeth
 
 
```
         set -euo pipefail
 
         if [ -z "~{numReads}" ]; then
             totalLines=$(zcat ~{fastqR1} | wc -l)
         else totalLines=$((~{numReads}*4))
         fi
         
         python3 -c "from math import ceil; print (int(ceil(($totalLines/4.0)/~{numChunk})*4))"
```
```
         set -euo pipefail
         slicer -i ~{fastqR} -l ~{chunkSize} --gzip 
```
```
         set -euo pipefail
         fastp \
             --stdout --thread ~{threads} \
             ~{fastpQ} ~{fastpq} ~{fastpu} ~{fastpn} ~{fastpL} ~{fastpl} ~{fastpA} ~{fastpG} \
             -i ~{read1} -I ~{read2} \
         | bwameth.py -p --threads ~{threads} --read-group ~{bwaReadGroup} --reference ~{bwaIndex} /dev/stdin \
         | samtools sort -o output.bam -@ ~{threads} -
```
```
         set -euo pipefail
 
         export JAVA_OPTS="-Xmx$(echo "scale=0; ~{jobMemory} * 0.8 / 1" | bc)G"
         java -jar ${PICARD_ROOT}/picard.jar \
         MergeSamFiles \
         I=~{sep=" I=" bams} \
         O=mergedChunks.bam \
         USE_THREADING=true \
         SORT_ORDER=coordinate
```
```
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
```

 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
