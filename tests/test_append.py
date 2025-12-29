# The idea here is to create two vcz datasets with different samples then we append on onto the other.
# For the moment don't worry about atomicity or transactions, as we can use icechunk later for that.
# Also, don't think about allele harmonisation yet.

# Create the VCF files, one with samples NA00001 and NA00002 and the other with NA00003

# bcftools view -s NA00001,NA00002 -O z tests/data/vcf/sample.vcf.gz > tests/data/vcf/sample-part1.vcf.gz
# bcftools view -s NA00003 -O z tests/data/vcf/sample.vcf.gz > tests/data/vcf/sample-part2.vcf.gz
# bcftools index -c tests/data/vcf/sample-part1.vcf.gz
# bcftools index -c tests/data/vcf/sample-part2.vcf.gz

# Create vcz datasets

# rm -rf ~/workspace/vczlib-poc/sample-part1.vcf.vcz; python -m bio2zarr vcf2zarr convert --variants-chunk-size 10 --samples-chunk-size 2 --force tests/data/vcf/sample-part1.vcf.gz ~/workspace/vczlib-poc/sample-part1.vcf.vcz
# rm -rf ~/workspace/vczlib-poc/sample-part2.vcf.vcz; python -m bio2zarr vcf2zarr convert --variants-chunk-size 10 --samples-chunk-size 2 --force tests/data/vcf/sample-part2.vcf.gz ~/workspace/vczlib-poc/sample-part2.vcf.vcz

# python -m vcztools view -H ~/workspace/vczlib-poc/sample-part1.vcf.vcz

# bcftools query -f '[%CHROM %POS %SAMPLE %GT\n]' ~/workspace/vcztools/tests/data/vcf/sample.vcf.gz 
# python -m vcztools query -f '[%CHROM %POS %SAMPLE %GT\n]' ~/workspace/vczlib-poc/sample-part1.vcf.vcz

# conda activate vczlib-poc-zarr-v2
# pytest -vs tests/test_append.py

# number of samples
# python -m vcztools query -l ~/workspace/vczlib-poc/sample-part1.vcf.vcz

# python -m vcztools view -H ~/workspace/vczlib-poc/sample-part1.vcf.vcz
# bcftools view -H tests/data/vcf/sample.vcf.gz

# check GT field
# diff <(python -m vcztools query -f '[%CHROM %POS %SAMPLE %GT\n]' ~/workspace/vczlib-poc/sample-part1.vcf.vcz) <(bcftools query -f '[%CHROM %POS %SAMPLE %GT\n]' tests/data/vcf/sample.vcf.gz)

import zarr

from vczlib import append

def test_append():
    append("sample-part1.vcf.vcz", "sample-part2.vcf.vcz")

def dims(arr):
    return arr.attrs["_ARRAY_DIMENSIONS"]

def test_find_samples_arrays():
    vcz = "sample-part1.vcf.vcz"

    root = zarr.open(vcz, mode="r")
    for var in root.keys():
        arr = root[var]
        if "samples" in dims(arr):
            print(var, dims(arr))
