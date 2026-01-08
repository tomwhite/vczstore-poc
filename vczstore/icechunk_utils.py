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
