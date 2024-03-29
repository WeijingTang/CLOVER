# CLOVER

This pipeline process two samples(one sample and it's corresponding Day0 sample) at one time, output an entropy ranking list(with efficiency) to each design.

# Parameters
sample(*fastq.gz), Day0_sample(*fastq.gz), sample_barcode_list(column one: line number, column two: barcode list), Day0_sample_barcode_list(column one: line number, column two: barcode list), whitelist, downsample_number

Please use zgrep -onP "T{6}[ATCG]{14}TTT[ATCG]" sample > sample.bcgrep.txt to generate sample_barcode_list

# Procedure
There are six steps in data_processing.py 
1.  Sample demutiplexing:
   it read the sample, check if the barcode in each squence in whitelist and put all sequences with the same whitelist into one fastq file
2.  Grep target region:
   to each sequence, keep barcode and target region.
3.  Sequence alignment:
   data_processing.py will call needleall_alignment.sh automatically and output needleall alignment results.
4.  Day0 filter:
   read all sample alignment results and remove the reads have the same Cigar String within the same design in Day0 alignment.
5.  Downsample:
   you can set a number, it will automatically downsample to that number.
6.  Entropy calculation:
   for each sample, there will be an entropy ranking list.

# Output
1. alignment results
2. entropy ranking list

Tip: please set you own path
