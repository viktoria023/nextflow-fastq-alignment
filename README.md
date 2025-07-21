# nextflow-fastq-alignment
Repository for aligning specified FASTQ reads from M. tuberculosis samples to a reference genome. Supported are long reads (ONT-Nanopore) and paired short reads (Illumina). Requires Nextflow and Conda to be installed.

**Installation**
```
git clone https://github.com/viktoria023/nextflow-fastq-alignment.git
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
    [--outputDir ./results]

# Long reads (Nanopore)
nextflow run read_alignment_pipeline.nf \
    --readType long \
    --inputFastq ./fastq_files \
    --referenceFasta ./reference.fasta \
    [--outputPrefix 'sample1'] \
    [--outputDir ./results]
```

**Input file structure**

Short reads
```
fastq_files/
├── sample1_1.fastq.gz
└── sample1_2.fastq.gz
```

Long reads
```
fastq_files/
└── sample1.fastq.gz
```

**Pipeline structure**
| Step  | Tool for short reads | Tool for long reads |
| ------------- | ------------- | ------------- |
| Alignment  | minimap2 `-ax sr`  | minimap2 `-ax map-ont` |
| Variant Calling  | FreeBayes  | longshot |
| Consensus  | bcftools consensus  | bcftools consensus |

