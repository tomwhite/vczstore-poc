import cubed
from cubed.array.update import append as cubed_append
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
    # Zarr format v2 has _ARRAY_DIMENSIONS, v3 has dedicated metadata
    return arr.attrs.get("_ARRAY_DIMENSIONS", None) or arr.metadata.dimension_names

INT_MISSING, INT_FILL = -1, -2

FLOAT32_MISSING, FLOAT32_FILL = np.array([0x7F800001, 0x7F800002], dtype=np.int32).view(
    np.float32
)
FLOAT32_MISSING_AS_INT32, FLOAT32_FILL_AS_INT32 = np.array(
    [0x7F800001, 0x7F800002], dtype=np.int32
)

# end from vcztools

def append(vcz1, vcz2):
    """Append vcz2 to vcz1 in place"""
    root1 = zarr.open(vcz1, mode="r+")

    cubed_arrays = []

    # append samples
    sample_id1 = root1["sample_id"]
    sample_id2 = cubed.from_zarr(vcz2, path="sample_id")
    c = cubed_append(sample_id1, sample_id2, axis=0)
    cubed_arrays.append(c)

    # append genotype fields
    for var in root1.keys():
        if var.startswith("call_"):
            arr = root1[var]
            arr2 = cubed.from_zarr(vcz2, path=var)
            c = cubed_append(arr, arr2, axis=1)
            cubed_arrays.append(c)

    # compute all arrays
    cubed.compute(*cubed_arrays, _return_in_memory_array=False)

    # consolidate metadata
    zarr.consolidate_metadata(vcz1)
