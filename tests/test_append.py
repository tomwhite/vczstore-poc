# Create the VCF files, one with samples NA00001 and NA00002 and the other with NA00003

# bcftools view -s NA00001,NA00002 --no-update -O z tests/data/vcf/sample.vcf.gz \
#  > tests/data/vcf/sample-part1.vcf.gz
# bcftools view -s NA00003 --no-update -O z tests/data/vcf/sample.vcf.gz \
#  > tests/data/vcf/sample-part2.vcf.gz
# bcftools index -c tests/data/vcf/sample-part1.vcf.gz
# bcftools index -c tests/data/vcf/sample-part2.vcf.gz

import pytest

from vczstore.zarr_impl import append as zarr_impl_append

try:
    from vczstore.cubed_impl import append as cubed_impl_append
except ImportError:
    cubed_impl_append = None

try:
    from vczstore.xarray_impl import append as xarray_impl_append
except ImportError:
    xarray_impl_append = None

from .utils import (
    compare_vcf_and_vcz,
    convert_vcf_to_vcz,
    convert_vcf_to_vcz_icechunk,
    run_vcztools,
)


@pytest.mark.parametrize(
    "append_function",
    [
        pytest.param(zarr_impl_append, id="zarr_impl"),
        pytest.param(
            cubed_impl_append,
            id="cubed_impl",
            marks=pytest.mark.skipif(
                cubed_impl_append is None, reason="requires cubed"
            ),
        ),
        pytest.param(
            xarray_impl_append,
            id="xarray_impl",
            marks=pytest.mark.skipif(
                xarray_impl_append is None, reason="requires xarray"
            ),
        ),
    ],
)
def test_append(tmp_path, append_function):
    print(tmp_path)

    vcz1 = convert_vcf_to_vcz("sample-part1.vcf.gz", tmp_path)
    vcz2 = convert_vcf_to_vcz("sample-part2.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1}")
    assert vcztools_out.strip() == "NA00001\nNA00002"

    append_function(vcz1, vcz2)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1}")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    # check equivalence with original VCF
    compare_vcf_and_vcz(
        tmp_path, "view --no-version", "sample.vcf.gz", "view --no-version", vcz1
    )


def test_append_icechunk(tmp_path):
    pytest.importorskip("icechunk")
    from icechunk import Repository

    from vczstore.icechunk_utils import make_icechunk_storage

    print(tmp_path)

    # note that vcz1 is in icechunk, but the dataset being appended, vcz2, needn't be
    vcz1 = convert_vcf_to_vcz_icechunk("sample-part1.vcf.gz", tmp_path)
    vcz2 = convert_vcf_to_vcz("sample-part2.vcf.gz", tmp_path)

    print(vcz1)
    print(vcz2)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1} --zarr-backend-storage icechunk")
    assert vcztools_out.strip() == "NA00001\nNA00002"

    icechunk_storage = make_icechunk_storage(vcz1)
    repo = Repository.open(icechunk_storage)

    with repo.transaction("main", message="append") as store:
        zarr_impl_append(store, vcz2)

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


def test_append_icechunk_cubed(tmp_path):
    pytest.importorskip("icechunk")
    pytest.importorskip("cubed")
    from icechunk import Repository

    from vczstore.cubed_impl import append as cubed_impl_append
    from vczstore.icechunk_utils import make_icechunk_storage

    print(tmp_path)

    # note that vcz1 is in icechunk, but the dataset being appended, vcz2, needn't be
    vcz1 = convert_vcf_to_vcz_icechunk("sample-part1.vcf.gz", tmp_path)
    vcz2 = convert_vcf_to_vcz("sample-part2.vcf.gz", tmp_path)

    print(vcz1)
    print(vcz2)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz1} --zarr-backend-storage icechunk")
    assert vcztools_out.strip() == "NA00001\nNA00002"

    icechunk_storage = make_icechunk_storage(vcz1)
    repo = Repository.open(icechunk_storage)

    session = repo.writable_session("main")
    fork = session.fork()
    store = fork.store

    from cubed import config

    with config.set({"spec.executor_name": "processes"}):
        merged_session = cubed_impl_append(store, vcz2, icechunk=True)

    session.merge(merged_session)
    session.commit("append")

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
