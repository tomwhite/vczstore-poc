import pytest

from vczstore import remove

from .utils import (
    check_removed_sample,
    compare_vcf_and_vcz,
    convert_vcf_to_vcz,
    convert_vcf_to_vcz_icechunk,
    run_vcztools,
)


def test_remove(tmp_path):
    print(tmp_path)

    vcz = convert_vcf_to_vcz("sample.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    remove(vcz, "NA00002")

    # TODO: note following requires the sample-mask branch of vcztools

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00003"

    # check equivalence with original VCF (with sample subsetting)
    compare_vcf_and_vcz(
        tmp_path,
        "view --no-version -s NA00001,NA00003",
        "sample.vcf.gz",
        "view --no-version",
        vcz,
    )


def test_remove_cubed(tmp_path):
    pytest.importorskip("cubed")

    print(tmp_path)

    vcz = convert_vcf_to_vcz("sample.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    from vczstore.cubed_impl import remove as remove_2

    remove_2(vcz, "NA00002")

    # TODO: note following requires the sample-mask branch of vcztools

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00003"

    # check equivalence with original VCF (with sample subsetting)
    compare_vcf_and_vcz(
        tmp_path,
        "view --no-version -s NA00001,NA00003",
        "sample.vcf.gz",
        "view --no-version",
        vcz,
    )


def test_remove_xarray(tmp_path):
    pytest.importorskip("xarray")
    from vczstore.xarray_impl import remove

    print(tmp_path)

    vcz = convert_vcf_to_vcz("sample.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    remove(vcz, "NA00002")

    # TODO: note following requires the sample-mask branch of vcztools

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz}")
    assert vcztools_out.strip() == "NA00001\nNA00003"

    # check equivalence with original VCF (with sample subsetting)
    compare_vcf_and_vcz(
        tmp_path,
        "view --no-version -s NA00001,NA00003",
        "sample.vcf.gz",
        "view --no-version",
        vcz,
    )

    # check sample values are missing
    check_removed_sample(vcz, "NA00002")


def test_remove_icechunk(tmp_path):
    pytest.importorskip("icechunk")
    from icechunk import Repository, Storage

    from vczstore.icechunk_utils import delete_previous_snapshots

    print(tmp_path)

    vcz = convert_vcf_to_vcz_icechunk("sample.vcf.gz", tmp_path)

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz} --zarr-backend-storage icechunk")
    assert vcztools_out.strip() == "NA00001\nNA00002\nNA00003"

    icechunk_storage = Storage.new_local_filesystem(str(vcz))
    repo = Repository.open(icechunk_storage)

    snapshots = [snapshot for snapshot in repo.ancestry(branch="main")]
    assert len(snapshots) == 2
    assert snapshots[0].message == "commit 1"
    assert snapshots[1].message == "Repository initialized"

    with repo.transaction("main", message="append") as store:
        remove(store, "NA00002", consolidate_metadata=False)

    delete_previous_snapshots(repo)

    snapshots = [snapshot for snapshot in repo.ancestry(branch="main")]
    assert len(snapshots) == 2
    # note that 'commit 1' has been deleted
    assert snapshots[0].message == "append"
    assert snapshots[1].message == "Repository initialized"

    # check samples query
    vcztools_out, _ = run_vcztools(f"query -l {vcz} --zarr-backend-storage icechunk")
    assert vcztools_out.strip() == "NA00001\nNA00003"

    # check equivalence with original VCF (with sample subsetting)
    compare_vcf_and_vcz(
        tmp_path,
        "view --no-version -s NA00001,NA00003",
        "sample.vcf.gz",
        "view --no-version --zarr-backend-storage icechunk",
        vcz,
    )
