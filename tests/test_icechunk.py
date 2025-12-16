from typing import Iterable

import numpy as np
import icechunk
import pytest
import tempfile
import zarr

from icechunk import Repository, Storage


@pytest.fixture
def icechunk_storage(tmpdir) -> "Storage":
    return Storage.new_local_filesystem(str(tmpdir))


def create_icechunk(a, icechunk_storage, /, *, dtype=None, chunks=None):
    if not isinstance(getattr(a, "shape", None), Iterable):
        # ensure blocks are arrays
        a = np.asarray(a, dtype=dtype)
    if dtype is None:
        dtype = a.dtype

    repo = Repository.create(storage=icechunk_storage)
    session = repo.writable_session("main")
    store = session.store

    group = zarr.group(store=store, overwrite=True)
    arr = group.create_array("a", shape=a.shape, dtype=dtype, chunks=chunks)

    arr[...] = a

    snapshot = session.commit("commit 1")

    return snapshot

def test_create_icechunk(icechunk_storage):
    snapshot = create_icechunk(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        icechunk_storage,
        chunks=(2, 2),
    )
    print(snapshot)

    # reopen store and check contents of array
    repo = Repository.open(icechunk_storage)

    print([ancestor for ancestor in repo.ancestry(branch="main")])

    session = repo.readonly_session(branch="main")
    store = session.store

    group = zarr.open_group(store=store, mode="r")
    np.testing.assert_array_equal(
        group["a"][:], np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    )

def test_update_icechunk(icechunk_storage):
    repo = Repository.create(icechunk_storage,
                             config=icechunk.RepositoryConfig(inline_chunk_threshold_bytes=0))

    a = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    with repo.transaction("main", message="create a") as store:
        root = zarr.open_group(store)
        arr = root.create_array("a", shape=a.shape, dtype=a.dtype, chunks=(2, 2))
        arr[...] = a

    with repo.transaction("main", message="update a") as store:
        root = zarr.open_group(store)
        root["a"][:, 1] = -1  # "delete" sample at index 1

    # reopen store and check contents of array
    repo = Repository.open(icechunk_storage)

    print("snapshots...")
    for snapshot in repo.ancestry(branch="main"):
        print(snapshot)

    session = repo.readonly_session(branch="main")
    store = session.store

    group = zarr.open_group(store=store, mode="r")
    np.testing.assert_array_equal(
        group["a"][:], np.array([[1, -1, 3], [4, -1, 6], [7, -1, 9]])
    )

    # expire old snapshots (see https://icechunk.io/en/stable/expiration/)
    current_snapshot = list(repo.ancestry(branch="main"))[0]
    expiry_time = current_snapshot.written_at
    expired = repo.expire_snapshots(older_than=expiry_time)
    print("expired...")
    print(expired)

    print("snapshots...")
    for snapshot in repo.ancestry(branch="main"):
        print(snapshot)

    session = repo.readonly_session(branch="main")
    store = session.store

    group = zarr.open_group(store=store, mode="r")
    np.testing.assert_array_equal(
        group["a"][:], np.array([[1, -1, 3], [4, -1, 6], [7, -1, 9]])
    )

    results = repo.garbage_collect(expiry_time)
    print("deleted...")
    print(results)


def test_append_icechunk(icechunk_storage):
    repo = Repository.create(icechunk_storage)

    a = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    with repo.transaction("main", message="create a") as store:
        root = zarr.open_group(store)
        arr = root.create_array("a", shape=a.shape, dtype=a.dtype, chunks=(2, 2))
        arr[...] = a

    b = -np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    with repo.transaction("main", message="append b") as store:
        root = zarr.open_group(store)
        arr = root["a"]
        arr.resize((3, 6))
        arr[:, 3:6] = b

    # reopen store and check contents of array
    repo = Repository.open(icechunk_storage)

    for snapshot in repo.ancestry(branch="main"):
        print(snapshot)

    session = repo.readonly_session(branch="main")
    store = session.store

    group = zarr.open_group(store=store, mode="r")
    np.testing.assert_array_equal(
        group["a"][:], np.concat((a, b), axis=1)
    )

# TODO: tag/release mechanism