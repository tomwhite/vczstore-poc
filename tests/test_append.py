# Create the VCF files, one with samples NA00001 and NA00002 and the other with NA00003

# bcftools view -s NA00001,NA00002 --no-update -O z tests/data/vcf/sample.vcf.gz \
#  > tests/data/vcf/sample-part1.vcf.gz
# bcftools view -s NA00003 --no-update -O z tests/data/vcf/sample.vcf.gz \
#  > tests/data/vcf/sample-part2.vcf.gz
# bcftools index -c tests/data/vcf/sample-part1.vcf.gz
# bcftools index -c tests/data/vcf/sample-part2.vcf.gz

import pytest

from vczstore import append

from .utils import (
    compare_vcf_and_vcz,
    convert_vcf_to_vcz,
    convert_vcf_to_vcz_icechunk,
    run_vcztools,
)


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
    compare_vcf_and_vcz(
        tmp_path, "view --no-version", "sample.vcf.gz", "view --no-version", vcz1
    )


def test_append_cubed(tmp_path):
    pytest.importorskip("cubed")

    print(tmp_path)

    vcz1 = convert_vcf_to_vcz("sample-part1.vcf.gz", tmp_path)
    vcz2 = convert_vcf_to_vcz("sample-part2.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1}")
    assert vcztools_out.strip() == "NA00001\nNA00002"

    from vczstore.cubed_impl import append as append_2

    append_2(vcz1, vcz2)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1}")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    # check equivalence with original VCF
    compare_vcf_and_vcz(
        tmp_path, "view --no-version", "sample.vcf.gz", "view --no-version", vcz1
    )


def test_append_icechunk(tmp_path):
    pytest.importorskip("icechunk")
    from icechunk import Repository, Storage

    print(tmp_path)

    # note that vcz1 is in icechunk, but the dataset being appended, vcz2, needn't be
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
    compare_vcf_and_vcz(
        tmp_path,
        "view --no-version",
        "sample.vcf.gz",
        "view --no-version --zarr-backend-storage icechunk",
        vcz1,
    )
