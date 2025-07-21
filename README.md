# nextflow-fastq-alignment
Repository for aligning specified FASTQ reads from M. tuberculosis samples to a reference genome. Supported are long reads (ONT-Nanopore) and paired short reads (Illumina).

**Installation**
```
git clone <repository-url>
cd nextflow-fastq-alignment
```

**Usage**
```
# Short reads (Illumina paired-end)
nextflow run read_alignment_pipeline.nf \
    --readType short \
    --inputFastq ./fastq_files \
    --referenceFasta ./reference.fasta \
    [--outputPrefix 'sample1'] \
    [--outputDir .results]

# Long reads (Nanopore)
nextflow run read_alignment_pipeline.nf \
    --readType long \
    --inputFastq ./fastq_files \
    --referenceFasta ./reference.fasta \
    [--outputPrefix 'sample1'] \
    [--outputDir .results]
```

**Pipeline structure**
| Step  | Tool for short reads | Tool for long reads |
| ------------- | ------------- |
| Alignment  | minimap2 `-ax sr`  | minimap2 `-ax map-ont` |
| Variant Calling  | FreeBayes  | longshot |
| Consensus  | bcftools consensus  | bcftools consensus |

