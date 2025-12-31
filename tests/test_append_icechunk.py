import pytest

from .utils import compare_vcf_and_vcz, convert_vcf_to_vcz, convert_vcf_to_vcz_icechunk, run_vcztools

from vczlib import append

icechunk = pytest.importorskip("icechunk")
from icechunk import Repository, Storage

def test_append(tmp_path):
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
        append(store, vcz2)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1} --zarr-backend-storage icechunk")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    # check equivalence with original VCF
    compare_vcf_and_vcz(tmp_path, "view --no-version", "sample.vcf.gz", "view --no-version --zarr-backend-storage icechunk", vcz1)
