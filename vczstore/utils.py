from vcztools.constants import FLOAT32_MISSING, INT_MISSING


def missing_val(arr):
    if arr.dtype.kind == "i":
        return INT_MISSING
    elif arr.dtype.kind == "f":
        return FLOAT32_MISSING
    elif arr.dtype.kind == "b":
        return False
    else:
        raise ValueError(f"unrecognised dtype: {arr.dtype}")
