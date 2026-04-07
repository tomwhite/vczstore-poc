import dataclasses

import numpy as np
import zarr
from bio2zarr.core import JsonDataclass
from bio2zarr.vcz import VcfZarrPartition
from vcztools.utils import search
from vcztools.vcf_writer import dims

from vczstore.utils import missing_val


@dataclasses.dataclass
class AppendWorkSummary(JsonDataclass):
    num_partitions: int
    num_variants: int


@dataclasses.dataclass
class RemoveWorkSummary(JsonDataclass):
    num_partitions: int
    num_variants: int


def append_init(vcz1, vcz2, target_num_partitions) -> AppendWorkSummary:
    root1 = zarr.open(vcz1, mode="r+")
    root2 = zarr.open(vcz2, mode="r")

    # calculate actual number of partitions
    # TODO: check all variants dims are the same
    num_variants = root1["variant_position"].shape[0]
    variants_chunk_size = root1["variant_position"].chunks[0]
    partitions = VcfZarrPartition.generate_partitions(
        num_variants, variants_chunk_size, target_num_partitions
    )
    num_partitions = len(partitions)

    # check preconditions
    n_variants1 = root1["variant_contig"].shape[0]
    n_variants2 = root2["variant_contig"].shape[0]
    if n_variants1 != n_variants2:
        raise ValueError(
            "Stores being appended must have same number of variants. "
            f"First has {n_variants1}, second has {n_variants2}"
        )
    # note that variant field precondition checks are made in append_partition
    # if any fail then the data may be left in an inconsistent state

    # append samples
    sample_id1 = root1["sample_id"]
    sample_id2 = root2["sample_id"]

    old_num_samples = sample_id1.shape[0]
    new_num_samples = old_num_samples + sample_id2.shape[0]
    new_shape = (new_num_samples,)
    sample_id1.resize(new_shape)
    sample_id1[old_num_samples:new_num_samples] = sample_id2[:]

    # resize genotype fields
    for var in root1.keys():
        if var.startswith("call_"):
            arr = root1[var]
            if arr.ndim == 2:
                new_shape = (arr.shape[0], new_num_samples)
                arr.resize(new_shape)
            elif arr.ndim == 3:
                new_shape = (arr.shape[0], new_num_samples, arr.shape[2])
                arr.resize(new_shape)
            else:
                raise ValueError("unsupported number of dims")

    # store append parameters in zarr attributes
    root1.attrs["append"] = {
        "num_partitions": num_partitions,
        "num_variants": num_variants,
        "variants_chunk_size": variants_chunk_size,
        "old_num_samples": old_num_samples,
        "new_num_samples": new_num_samples,
    }

    return AppendWorkSummary(num_partitions, num_variants)


def append_partition(vcz1, vcz2, partition_index):
    root1 = zarr.open(vcz1, mode="r+")
    root2 = zarr.open(vcz2, mode="r")

    append_attrs = root1.attrs["append"]
    num_partitions = int(append_attrs["num_partitions"])
    num_variants = int(append_attrs["num_variants"])
    variants_chunk_size = int(append_attrs["variants_chunk_size"])
    old_num_samples = int(append_attrs["old_num_samples"])
    new_num_samples = int(append_attrs["new_num_samples"])

    partitions = VcfZarrPartition.generate_partitions(
        num_variants, variants_chunk_size, num_partitions
    )
    partition = partitions[partition_index]
    partition_sel = slice(partition.start, partition.stop)

    # check preconditions
    for field in ("contig_id", "variant_contig", "variant_position", "variant_allele"):
        values1 = root1[field][partition_sel]
        values2 = root2[field][partition_sel]
        if np.any(values1 != values2):
            raise ValueError(
                f"Stores being appended must have same values for field '{field}'"
            )

    # append genotype fields
    for var in root1.keys():
        if var.startswith("call_"):
            arr = root1[var]
            sl = slice(partition.start, partition.stop)
            if arr.ndim == 2:
                arr[partition_sel, old_num_samples:new_num_samples] = root2[var][
                    sl, ...
                ]
            elif arr.ndim == 3:
                arr[partition_sel, old_num_samples:new_num_samples, :] = root2[var][
                    sl, ...
                ]
            else:
                raise ValueError("unsupported number of dims")


def append_finalise(vcz1, vcz2):
    root = zarr.open(vcz1, mode="r+")
    del root.attrs["append"]


def remove_init(vcz, sample_id, target_num_partitions) -> RemoveWorkSummary:
    root = zarr.open(vcz, mode="r+")

    # calculate actual number of partitions
    # TODO: check all variants dims are the same
    num_variants = root["variant_position"].shape[0]
    variants_chunk_size = root["variant_position"].chunks[0]
    partitions = VcfZarrPartition.generate_partitions(
        num_variants, variants_chunk_size, target_num_partitions
    )
    num_partitions = len(partitions)

    all_samples = root["sample_id"][:]

    # find index of sample to remove
    unknown_samples = np.setdiff1d(sample_id, all_samples)
    if len(unknown_samples) > 0:
        raise ValueError(f"unrecognised sample: {sample_id}")
    selection = search(all_samples, sample_id)

    # overwrite sample data
    root["sample_id"][selection] = ""

    # store remove parameters in zarr attributes
    root.attrs["remove"] = {
        "num_partitions": num_partitions,
        "num_variants": num_variants,
        "variants_chunk_size": variants_chunk_size,
        "sample_index": int(selection),
    }

    return RemoveWorkSummary(num_partitions, num_variants)


def remove_partition(vcz, partition_index):
    root = zarr.open(vcz, mode="r+")

    remove_attrs = root.attrs["remove"]
    num_partitions = int(remove_attrs["num_partitions"])
    num_variants = int(remove_attrs["num_variants"])
    variants_chunk_size = int(remove_attrs["variants_chunk_size"])
    sample_index = int(remove_attrs["sample_index"])

    partitions = VcfZarrPartition.generate_partitions(
        num_variants, variants_chunk_size, num_partitions
    )
    partition = partitions[partition_index]
    partition_sel = slice(partition.start, partition.stop)

    # overwrite call variables
    for var in root.keys():
        arr = root[var]
        if (
            var.startswith("call_")
            and dims(arr)[0] == "variants"
            and dims(arr)[1] == "samples"
        ):
            root[var][partition_sel, sample_index, ...] = missing_val(arr)


def remove_finalise(vcz):
    root = zarr.open(vcz, mode="r+")
    del root.attrs["remove"]
