import numpy as np
import zarr
from vcztools.utils import search
from vcztools.vcf_writer import dims

from vczstore.allele_harmonisation import harmonise_alleles, remap_gt
from vczstore.utils import missing_val


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
    sample_id1.resize(new_shape)
    sample_id1[old_num_samples:new_num_samples] = sample_id2[:]

    # append genotype fields

    for var in root1.keys():
        if var.startswith("call_"):
            arr = root1[var]
            if arr.ndim == 2:
                new_shape = (arr.shape[0], new_num_samples)
                arr.resize(new_shape)
                arr[:, old_num_samples:new_num_samples] = root2[var][:]
            elif arr.ndim == 3:
                new_shape = (arr.shape[0], new_num_samples, arr.shape[2])
                arr.resize(new_shape)
                arr[:, old_num_samples:new_num_samples, :] = root2[var][:]
            else:
                raise ValueError("unsupported number of dims")

    # consolidate metadata (if supported)
    try:
        zarr.consolidate_metadata(vcz1)
    except TypeError:
        # store doesn't support consolidated metadata, that's OK
        pass


def append_harmonise(vcz1, vcz2):
    """Append vcz2 to vcz1 in place, harmonising alleles"""
    root1 = zarr.open(vcz1, mode="r+")
    # TODO: we'll update vcz2 in place to fix alleles (not sure this is OK?)
    root2 = zarr.open(vcz2, mode="r+")

    # append samples
    sample_id1 = root1["sample_id"]
    sample_id2 = root2["sample_id"]

    old_num_samples = sample_id1.shape[0]
    new_num_samples = old_num_samples + sample_id2.shape[0]
    new_shape = (new_num_samples,)
    sample_id1.resize(new_shape)
    sample_id1[old_num_samples:new_num_samples] = sample_id2[:]

    # harmonise alleles
    variant_allele = root1["variant_allele"]
    variant_allele_new = root2["variant_allele"]
    variant_allele_updated, variant_allele_new_mapping = harmonise_alleles(
        variant_allele, variant_allele_new
    )
    variant_allele.resize(variant_allele_updated.shape)
    variant_allele[:] = variant_allele_updated

    # append genotype fields

    for var in root1.keys():
        if var == "call_genotype":
            arr = root1[var]
            new_shape = (arr.shape[0], new_num_samples, arr.shape[2])
            arr.resize(new_shape)

            # TODO: consider case where ploidy dim needs resizing (expanding)
            gt = root2[var][:]
            remap_gt(gt, variant_allele_new_mapping)

            arr[:, old_num_samples:new_num_samples, :] = gt
        elif var.startswith("call_"):
            arr = root1[var]
            if arr.ndim == 2:
                new_shape = (arr.shape[0], new_num_samples)
                arr.resize(new_shape)
                arr[:, old_num_samples:new_num_samples] = root2[var][:]
            elif arr.ndim == 3:
                new_shape = (arr.shape[0], new_num_samples, arr.shape[2])
                arr.resize(new_shape)
                arr[:, old_num_samples:new_num_samples, :] = root2[var][:]
            else:
                raise ValueError("unsupported number of dims")

    # consolidate metadata (if supported)
    try:
        zarr.consolidate_metadata(vcz1)
    except TypeError:
        # store doesn't support consolidated metadata, that's OK
        pass


def remove(vcz, sample_id):
    """Remove a sample from vcz and overwrite with missing data"""
    root = zarr.open(vcz, mode="r+")
    all_samples = root["sample_id"][:]

    # find index of sample to remove
    unknown_samples = np.setdiff1d(sample_id, all_samples)
    if len(unknown_samples) > 0:
        raise ValueError(f"unrecognised sample: {sample_id}")
    selection = search(all_samples, sample_id)

    # overwrite sample data
    root["sample_id"][selection] = ""
    for var in root.keys():
        arr = root[var]
        if (
            var.startswith("call_")
            and dims(arr)[0] == "variants"
            and dims(arr)[1] == "samples"
        ):
            root[var][:, selection, ...] = missing_val(arr)

    # consolidate metadata (if supported)
    try:
        zarr.consolidate_metadata(vcz)
    except TypeError:
        # store doesn't support consolidated metadata, that's OK
        pass
