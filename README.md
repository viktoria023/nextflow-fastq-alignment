# nextflow-fastq-alignment
Repository for aligning specified FASTQ reads from M. tuberculosis samples to a reference genome. Supported are long reads (ONT-Nanopore) and paired short reads (Illumina). Requires Nextflow and Conda to be installed.

**Installation**
```bash
git clone https://github.com/viktoria023/nextflow-fastq-alignment.git
cd nextflow-fastq-alignment
```

**Usage**
```bash
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

**Minimal testing setup**
A minimal integration testing script was added, allowing to test run the pipeline for both read types using minimal synthetic test data and reference (~16KB):

```bash
# Run automated tests for both read types
chmod +x test_pipeline.sh
./test_pipeline.sh
```

Expected output:
```
📊 TEST SUMMARY
===============
✅ Short reads: PASSED
✅ Long reads: PASSED
```