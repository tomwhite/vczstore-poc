# The idea here is to create a vcz dataset then remove a sample.
# For the moment don't worry about atomicity or transactions, as we can use icechunk later for that.
# But we do need to consider masking and recomputing fields like AC.

# NOTE: the below is now out of date as the test generates its own data

# conda activate vczlib-poc-zarr-v2
# rm -rf sample.vcf.vcz; cp -r ~/workspace/vcztools/vcz_test_cache/sample.vcf.vcz sample.vcf.vcz
# pytest -vs tests/test_remove.py

# then in vcztools dir, delete-mask branch, vcztools-3.12 venv
# python -m vcztools query -f '[%SAMPLE %GT %DP\n]' -s NA00001,NA00003 ~/workspace/vczlib-poc/sample.vcf.vcz
# python -m vcztools view ~/workspace/vczlib-poc/sample.vcf.vcz
# python -m vcztools view -s NA00001,NA00003 ~/workspace/vczlib-poc/sample.vcf.vcz

# number of samples
# python -m vcztools query -l ~/workspace/vczlib-poc/sample.vcf.vcz

import pytest
from .utils import compare_vcf_and_vcz, convert_vcf_to_vcz, run_vcztools
from vczlib import remove

def test_remove(tmp_path):
    print(tmp_path)

    vcz = convert_vcf_to_vcz("sample.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    remove(vcz, "NA00002")

    # TODO: note following requires the delete-mask branch of vcztools

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00003"

    # check equivalence with original VCF (with sample subsetting)
    compare_vcf_and_vcz(tmp_path, "view --no-version -s NA00001,NA00003", "sample.vcf.gz", "view --no-version", vcz)


def test_remove_cubed(tmp_path):
    pytest.importorskip("cubed")

    print(tmp_path)

    vcz = convert_vcf_to_vcz("sample.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    from vczlib.cubed_impl import remove as remove_2
    remove_2(vcz, "NA00002")

    # TODO: note following requires the delete-mask branch of vcztools

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00003"

    # check equivalence with original VCF (with sample subsetting)
    compare_vcf_and_vcz(tmp_path, "view --no-version -s NA00001,NA00003", "sample.vcf.gz", "view --no-version", vcz)