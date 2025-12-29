# The idea here is to create a vcz dataset then remove a sample.
# For the moment don't worry about atomicity or transactions, as we can use icechunk later for that.
# But we do need to consider masking and recomputing fields like AC.

# conda activate vczlib-poc-zarr-v2
# rm -rf sample.vcf.vcz; cp -r ~/workspace/vcztools/vcz_test_cache/sample.vcf.vcz sample.vcf.vcz
# pytest -vs tests/test_remove.py

# then in vcztools dir, delete-mask branch, vcztools-3.12 venv
# python -m vcztools query -f '[%SAMPLE %GT %DP\n]' -s NA00001,NA00003 ~/workspace/vczlib-poc/sample.vcf.vcz
# python -m vcztools view ~/workspace/vczlib-poc/sample.vcf.vcz
# python -m vcztools view -s NA00001,NA00003 ~/workspace/vczlib-poc/sample.vcf.vcz

# number of samples
# python -m vcztools query -l ~/workspace/vczlib-poc/sample.vcf.vcz

from vczlib import remove

def test_remove():
    remove("sample.vcf.vcz", "NA00002")