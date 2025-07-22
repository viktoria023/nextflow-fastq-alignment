#!/bin/bash

echo "Running MTB consensus pipeline test"

# Check test data exists
if [ ! -f "test_data/reference/H37Rv.fasta" ] || [ ! -f "test_data/fastq/short_read/test_short_1.fastq.gz" ] || [ ! -f "test_data/fastq/long_read/test_long.fastq.gz" ]; then
    echo "❌ ERROR: Test data not found!"
    echo "Expected files:"
    echo "  - test_data/reference/H37Rv.fasta"
    echo "  - test_data/fastq/short_read/test_short_1.fastq.gz"
    echo "  - test_data/fastq/short_read/test_short_2.fastq.gz"
    echo "  - test_data/fastq/long_read/test_long.fastq.gz"
    exit 1
fi

echo "✅ Test data found!"

# Run the pipeline for short reads
echo ""
echo "🚀 Test 1: Running pipeline test for short reads"
nextflow run read_alignment_pipeline.nf \
    --readType short \
    --inputFastq ./test_data/fastq/short_read \
    --referenceFasta ./test_data/reference/H37Rv.fasta \
    --outputPrefix test_short

# Check short reads results
if [ -f "results/consensus/test_short_test_short.consensus.fasta" ]; then
    echo "✅ SHORT reads test PASSED!"
    short_test_passed=true
else
    echo "❌ SHORT reads test FAILED!"
    short_test_passed=false
fi

# Clean up for next test
rm -rf results work .nextflow*

# Run the pipeline for long reads
echo ""
echo "🚀 Test 2: Running pipeline test for long reads"
nextflow run read_alignment_pipeline.nf \
    --readType long \
    --inputFastq ./test_data/fastq/long_read \
    --referenceFasta ./test_data/reference/H37Rv.fasta \
    --outputPrefix test_long

# Check long reads results
if [ -f "results/consensus/test_long_test_long.consensus.fasta" ]; then
    echo "✅ LONG reads test PASSED!"
    long_test_passed=true
else
    echo "❌ LONG reads test FAILED!"
    long_test_passed=false
fi

# Clean up working directory
rm -rf results work .nextflow*

# Final summary
echo ""
echo "📊 TEST SUMMARY"
echo "==============="
if [ "$short_test_passed" = true ]; then
    echo "✅ Short reads: PASSED"
else
    echo "❌ Short reads: FAILED"
fi

if [ "$long_test_passed" = true ]; then
    echo "✅ Long reads: PASSED"
else
    echo "❌ Long reads: FAILED"
fi