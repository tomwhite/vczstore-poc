# The idea here is to create two vcz datasets with different samples then we append on onto the other.
# For the moment don't worry about atomicity or transactions, as we can use icechunk later for that.
# Also, don't think about allele harmonisation.

# Create the VCF files, one with samples NA00001 and NA00002 and the other with NA00003

# bcftools view -s NA00001,NA00002 -O z tests/data/vcf/sample.vcf.gz > tests/data/vcf/sample-part1.vcf.gz
# bcftools view -s NA00003 -O z tests/data/vcf/sample.vcf.gz > tests/data/vcf/sample-part2.vcf.gz
# bcftools index -c tests/data/vcf/sample-part1.vcf.gz
# bcftools index -c tests/data/vcf/sample-part2.vcf.gz

# Create vcz datasets

# rm -rf ~/workspace/vczlib-poc/sample-part1.vcf.vcz; python -m bio2zarr vcf2zarr convert --variants-chunk-size 10 --samples-chunk-size 2 --force tests/data/vcf/sample-part1.vcf.gz ~/workspace/vczlib-poc/sample-part1.vcf.vcz
# rm -rf ~/workspace/vczlib-poc/sample-part2.vcf.vcz; python -m bio2zarr vcf2zarr convert --variants-chunk-size 10 --samples-chunk-size 2 --force tests/data/vcf/sample-part2.vcf.gz ~/workspace/vczlib-poc/sample-part2.vcf.vcz

# python -m vcztools view ~/workspace/vczlib-poc/sample-part1.vcf.vcz

# conda activate vczlib-poc-zarr-v2
# pytest -vs tests/test_append.py

# number of samples
# python -m vcztools query -l ~/workspace/vczlib-poc/sample-part1.vcf.vcz

# python -m vcztools view ~/workspace/vczlib-poc/sample-part1.vcf.vcz

import zarr

def append(vcz1, vcz2):
    """Append vcz2 to vcz1 in place"""
    root1 = zarr.open(vcz1, mode="r+")
    root2 = zarr.open(vcz2, mode="r")

    # append samples
    sample_id1 = root1["sample_id"]
    sample_id2 = root2["sample_id"]

    old_num_samples = sample_id1.shape[0]
    new_num_samples = old_num_samples + sample_id2.shape[0]
    new_shape = (new_num_samples,)
    print(new_shape)
    sample_id1.resize(new_shape)
    sample_id1[old_num_samples:new_num_samples] = sample_id2[:]


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
