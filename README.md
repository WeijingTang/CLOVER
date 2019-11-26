# CLOVER

This pipeline process two samples(one sample and it's corresponding Day0 sample) at one time, output an entropy ranking list(with efficiency) to each design.

# Parameters
sample(*fastq.gz), Day0_sample(*fastq.gz), sample_barcode_list(column one: line number, column two: barcode list), Day0_sample_barcode_list(column one: line number, column two: barcode list), whitelist, downsample_number

# Procedure
there are four steps in data_processing.py 
1. sample demutiplexing
   it read the sample, check if the barcode in each squence in whitelist and put all sequences with the same whitelist into one fastq file
2. grep target region
   to each sequence, keep barcode and target region.
3. sequence alignment
   data_processing.py will call needleall_alignment.sh automatically and output needleall alignment result.
4. Day0 filter
   read all sample alignment results and remove the reads have the same Cigar String within the same design.
5. downsample
   you can set a number, it will automatically downsample to that number.
6. entropy calculation
   for each sample, there will be an entropy ranking list.

# Output
1. alignment results
2. entropy ranking list

tip: please set you own path
