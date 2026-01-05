import cubed
import numpy as np
import zarr
from cubed.array.update import append as cubed_append
from cubed.array.update import set_scalar as cubed_set
from vcztools.constants import FLOAT32_MISSING, INT_MISSING
from vcztools.utils import search
from vcztools.vcf_writer import dims


def append(vcz1, vcz2):
    """Append vcz2 to vcz1 in place"""
    root1 = zarr.open(vcz1, mode="r+")

    cubed_arrays = []

    # append samples
    sample_id1 = cubed.from_zarr(vcz1, path="sample_id", mode="r+")
    sample_id2 = cubed.from_zarr(vcz2, path="sample_id")
    c = cubed_append(sample_id1, sample_id2, axis=0)
    cubed_arrays.append(c)

    # append genotype fields
    for var in root1.keys():
        if var.startswith("call_"):
            arr = cubed.from_zarr(vcz1, path=var, mode="r+")
            arr2 = cubed.from_zarr(vcz2, path=var)
            c = cubed_append(arr, arr2, axis=1)
            cubed_arrays.append(c)

    # compute all arrays
    cubed.compute(*cubed_arrays, _return_in_memory_array=False)

    # consolidate metadata
    zarr.consolidate_metadata(vcz1)


def missing_val(arr):
    if arr.dtype.kind == "i":
        return INT_MISSING
    elif arr.dtype.kind == "f":
        return FLOAT32_MISSING
    elif arr.dtype.kind == "b":
        return False
    else:
        raise ValueError(f"unrecognised dtype: {arr.dtype}")


def remove(vcz, sample_id):
    root = zarr.open(vcz, mode="r+")
    all_samples = root["sample_id"][:]

    # find index of sample to remove
    unknown_samples = np.setdiff1d(sample_id, all_samples)
    if len(unknown_samples) > 0:
        raise ValueError(f"unrecognised sample: {sample_id}")
    selection = search(all_samples, sample_id)
    sample_id_mask = np.zeros(all_samples.shape, dtype=bool)
    sample_id_mask[selection] = True

    # create or update the sample mask
    # TODO: the mask should be a part of bio2zarr eventually
    if "sample_id_mask" not in root:
        dimension_names = ["samples"]
        array = root.array(
            "sample_id_mask",
            data=sample_id_mask,
            shape=sample_id_mask.shape,
            chunks=sample_id_mask.shape,
            dtype=sample_id_mask.dtype,
            # TODO: compressor or codecs?
            # TODO: dimension_names for v3
        )
        array.attrs["_ARRAY_DIMENSIONS"] = dimension_names
    else:
        root["sample_id_mask"] |= sample_id_mask

    # overwrite sample data
    cubed_arrays = []
    for var in root.keys():
        arr = root[var]
        if (
            var.startswith("call_")
            and dims(arr)[0] == "variants"
            and dims(arr)[1] == "samples"
        ):
            cubed_arr = cubed.from_zarr(vcz, path=var, mode="r+")
            c = cubed_set(
                cubed_arr, (slice(None), selection, Ellipsis), missing_val(arr)
            )
            cubed_arrays.append(c)
            # root[var][:, selection, ...] = missing_val(arr)

    # compute all arrays
    cubed.compute(*cubed_arrays, _return_in_memory_array=False)

    # TODO: recalculate variant_AC, variant_AN
    # see _compute_info_fields in vcztools

    # consolidate metadata (may not be needed if sample_id_mask was already present)
    zarr.consolidate_metadata(vcz)
