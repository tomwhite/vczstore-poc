# The idea here is to create a vcz dataset then remove a sample.
# For the moment don't worry about atomicity or transactions, as we can use icechunk later for that.
# But we do need to consider masking and recomputing fields like AC.

# rm -rf sample.vcf.vcz; cp -r ~/workspace/vcztools/vcz_test_cache/sample.vcf.vcz sample.vcf.vcz

# python -m vcztools query -f '[%SAMPLE %GT %DP\n]' -s NA00001,NA00003 ~/workspace/vczlib-poc/sample.vcf.vcz
# python -m vcztools view ~/workspace/vczlib-poc/sample.vcf.vcz
# python -m vcztools view -s NA00001,NA00003 ~/workspace/vczlib-poc/sample.vcf.vcz

# number of samples
# python -m vcztools query -l ~/workspace/vczlib-poc/sample.vcf.vcz

import numpy as np
import zarr

# from vcztools
def search(a, v):
    """
    Finds the indices into an array a corresponding to the elements in v.
    The behaviour is undefined if any elements in v are not in a.
    """
    sorter = np.argsort(a)
    rank = np.searchsorted(a, v, sorter=sorter)
    return sorter[rank]

def dims(arr):
    return arr.attrs["_ARRAY_DIMENSIONS"]

INT_MISSING, INT_FILL = -1, -2

FLOAT32_MISSING, FLOAT32_FILL = np.array([0x7F800001, 0x7F800002], dtype=np.int32).view(
    np.float32
)
FLOAT32_MISSING_AS_INT32, FLOAT32_FILL_AS_INT32 = np.array(
    [0x7F800001, 0x7F800002], dtype=np.int32
)

# end from vcztools

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
    if "sample_id_delete" not in root:
        root.array(
            "sample_id_delete",
            data=sample_id_delete,
            shape=sample_id_delete.shape,
            chunks=sample_id_delete.shape,
            dtype=sample_id_delete.dtype,
            # TODO: compressor?
        )
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

def test_remove():
    remove("sample.vcf.vcz", "NA00002")