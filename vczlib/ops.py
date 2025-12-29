import zarr

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
                raise ValueError("unsupport number of dims")
