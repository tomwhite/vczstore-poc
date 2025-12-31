import numpy as np
import zarr

from vcztools.constants import INT_MISSING, FLOAT32_MISSING
from vcztools.utils import search
from vcztools.vcf_writer import dims

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

    # consolidate metadata
    # TODO: icechunk doesn't support consolidated metadata
    # zarr.consolidate_metadata(vcz1)

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
    sample_id_delete = np.zeros(all_samples.shape, dtype=bool)
    sample_id_delete[selection] = True

    # create or update the delete mask
    # TODO: the mask should be a part of bio2zarr eventually
    if "sample_id_delete" not in root:
        dimension_names = ["samples"]
        array = root.array(
            "sample_id_delete",
            data=sample_id_delete,
            shape=sample_id_delete.shape,
            chunks=sample_id_delete.shape,
            dtype=sample_id_delete.dtype,
            # TODO: compressor or codecs?
            # TODO: dimension_names for v3
        )
        array.attrs["_ARRAY_DIMENSIONS"] = dimension_names
    else:
        root["sample_id_delete"] |= sample_id_delete

    # overwrite sample data
    for var in root.keys():
        arr = root[var]
        if (
            var.startswith("call_")
            and dims(arr)[0] == "variants"
            and dims(arr)[1] == "samples"
        ):
            root[var][:, selection, ...] = missing_val(arr)

    # TODO: recalculate variant_AC, variant_AN
    # see _compute_info_fields in vcztools

    # consolidate metadata (may not be needed if sample_id_delete was already present)
    zarr.consolidate_metadata(vcz)