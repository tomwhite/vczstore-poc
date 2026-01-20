from pathlib import Path
from urllib.parse import urlparse

from zarr.core.sync import sync
from zarr.storage._common import make_store


def make_icechunk_storage(file_or_url):
    """Convert a file or URL to an Icechunk Storage object."""
    import icechunk as ic

    if isinstance(file_or_url, str):
        if "://" not in file_or_url:  # local path
            return ic.Storage.new_local_filesystem(file_or_url)
        elif file_or_url.startswith("s3://"):
            url_parsed = urlparse(file_or_url)
            return ic.s3_storage(
                bucket=url_parsed.netloc,
                prefix=url_parsed.path.lstrip("/"),
                from_env=True,
            )
        else:
            raise ValueError(f"Unsupported URL for icechunk: {file_or_url}")
    elif isinstance(file_or_url, Path):
        path = file_or_url.resolve()  # make absolute
        return ic.Storage.new_local_filesystem(str(path))
    else:
        raise TypeError(f"Unsupported URL type for icechunk: {type(file_or_url)}")


def delete_previous_snapshots(repo, branch="main"):
    """
    Delete all previous snapshots except the current one
    to avoid retaining data that has been explicitly removed.
    """
    # see https://icechunk.io/en/stable/expiration/

    current_snapshot = list(repo.ancestry(branch=branch))[0]
    expiry_time = current_snapshot.written_at
    repo.expire_snapshots(older_than=expiry_time)
    repo.garbage_collect(expiry_time)


# inspired by commit f3c123d3a2a94b7f14bc995e3897ee6acc9acbd1 in zarr-python
def copy_store(source, dest):
    from zarr.core.buffer.core import default_buffer_prototype
    from zarr.testing.stateful import SyncStoreWrapper

    # ensure source and dest are both stores
    source = sync(make_store(source))
    dest = sync(make_store(dest))

    s = SyncStoreWrapper(source)
    d = SyncStoreWrapper(dest)
    # need reverse=True to create zarr.json before chunks (otherwise icechunk complains)
    for source_key in sorted(s.list(), reverse=True):
        buffer = s.get(source_key, default_buffer_prototype())
        d.set(source_key, buffer)
