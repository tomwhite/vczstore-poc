# The idea here is to create two vcz datasets with different samples then we append on onto the other.
# For the moment don't worry about atomicity or transactions, as we can use icechunk later for that.
# Also, don't think about allele harmonisation yet.

# Create the VCF files, one with samples NA00001 and NA00002 and the other with NA00003

# bcftools view -s NA00001,NA00002 --no-update -O z tests/data/vcf/sample.vcf.gz > tests/data/vcf/sample-part1.vcf.gz
# bcftools view -s NA00003 --no-update -O z tests/data/vcf/sample.vcf.gz > tests/data/vcf/sample-part2.vcf.gz
# bcftools index -c tests/data/vcf/sample-part1.vcf.gz
# bcftools index -c tests/data/vcf/sample-part2.vcf.gz

# NOTE: the below is now out of date as the test generates its own data

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

import pytest

from .utils import compare_vcf_and_vcz, convert_vcf_to_vcz, convert_vcf_to_vcz_icechunk, run_vcztools
import zarr

from vczlib import append

def test_append(tmp_path):
    print(tmp_path)

    vcz1 = convert_vcf_to_vcz("sample-part1.vcf.gz", tmp_path)
    vcz2 = convert_vcf_to_vcz("sample-part2.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1}")
    assert vcztools_out.strip() == "NA00001\nNA00002"

    append(vcz1, vcz2)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1}")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    # check equivalence with original VCF
    compare_vcf_and_vcz(tmp_path, "view --no-version", "sample.vcf.gz", "view --no-version", vcz1)


def test_append_cubed(tmp_path):
    pytest.importorskip("cubed")

    print(tmp_path)

    vcz1 = convert_vcf_to_vcz("sample-part1.vcf.gz", tmp_path)
    vcz2 = convert_vcf_to_vcz("sample-part2.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1}")
    assert vcztools_out.strip() == "NA00001\nNA00002"

    from vczlib.cubed_impl import append as append_2
    append_2(vcz1, vcz2)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1}")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    # check equivalence with original VCF
    compare_vcf_and_vcz(tmp_path, "view --no-version", "sample.vcf.gz", "view --no-version", vcz1)


def test_append_icechunk(tmp_path):
    pytest.importorskip("icechunk")
    from icechunk import Repository, Storage

    print(tmp_path)

    # note that vcz1 is in icechunk, but the dataset being appended, vcz2, doesn't need to be
    vcz1 = convert_vcf_to_vcz_icechunk("sample-part1.vcf.gz", tmp_path)
    vcz2 = convert_vcf_to_vcz("sample-part2.vcf.gz", tmp_path)

    print(vcz1)
    print(vcz2)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1} --zarr-backend-storage icechunk")
    assert vcztools_out.strip() == "NA00001\nNA00002"

    icechunk_storage = Storage.new_local_filesystem(str(vcz1))
    repo = Repository.open(icechunk_storage)

    with repo.transaction("main", message="append") as store:
        append(store, vcz2, consolidate_metadata=False)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1} --zarr-backend-storage icechunk")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    # check equivalence with original VCF
    compare_vcf_and_vcz(tmp_path, "view --no-version", "sample.vcf.gz", "view --no-version --zarr-backend-storage icechunk", vcz1)