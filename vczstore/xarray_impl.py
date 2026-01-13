import numpy as np
import xarray as xr
from vcztools.constants import FLOAT32_MISSING, INT_MISSING
from vcztools.utils import search


def append(vcz1, vcz2):
    """Append vcz2 to vcz1 in place"""
    ds2 = xr.open_zarr(vcz2)
    ds2.to_zarr(vcz1, append_dim="samples")


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
    # find index of sample to remove
    ds = xr.open_zarr(vcz)
    all_samples = ds["sample_id"].values[:]
    unknown_samples = np.setdiff1d(sample_id, all_samples)
    if len(unknown_samples) > 0:
        raise ValueError(f"unrecognised sample: {sample_id}")
    index = search(all_samples, sample_id)
    sl = slice(index, index + 1)

    # update all the samples variables
    ds = xr.open_zarr(vcz)
    ds2 = ds.sel(samples=ds.sample_id.isin([sample_id]))
    non_sample_vars = [
        var
        for var in ds.data_vars
        if "samples" not in ds[var].dims or var == "sample_id"
    ]
    ds2 = ds2.drop_vars(non_sample_vars)
    for var in ds2.data_vars:
        da = ds2[var]
        da[...] = missing_val(da)
    region = {"samples": sl}
    ds2.to_zarr(vcz, region=region)

    # add sample_id_mask variable to dataset
    ds = xr.open_zarr(vcz)
    mask = np.zeros(ds.sizes["samples"], dtype=bool)
    mask[sl] = True
    sample_id_mask = xr.DataArray(mask, dims=("samples"))
    ds2 = xr.Dataset({"sample_id_mask": sample_id_mask})
    ds2.to_zarr(vcz, mode="a")
