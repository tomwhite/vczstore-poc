from contextlib import nullcontext

import click


class NaturalOrderGroup(click.Group):
    """
    List commands in the order they are provided in the help text.
    """

    def list_commands(self, ctx):
        return self.commands.keys()


impl = click.option(
    "-i",
    "--impl",
    type=str,
    default=None,
    help="Implementation to use; one of 'zarr' (default), 'cubed', or 'xarray'.",
)

zarr_backend_storage = click.option(
    "--zarr-backend-storage",
    type=str,
    default=None,
    help="Zarr backend storage to use; one of 'fsspec' (default) or 'icechunk'.",
)


@click.command()
@click.argument("vcz1", type=click.Path())
@click.argument("vcz2", type=click.Path())
@impl
@zarr_backend_storage
def append(vcz1, vcz2, impl, zarr_backend_storage):
    """Append vcz2 to vcz1 in place"""
    if impl is None or impl == "zarr":
        from vczstore.zarr_impl import append as append_function
    elif impl == "cubed":
        from vczstore.cubed_impl import append as append_function
    elif impl == "xarray":
        from vczstore.xarray_impl import append as append_function
    else:
        raise ValueError(f"Unrecognised impl: {impl}")

    if zarr_backend_storage == "icechunk":
        from vczstore.icechunk_utils import icechunk_transaction

        cm = icechunk_transaction(vcz1, "main", message="append")
    else:
        cm = nullcontext(vcz1)
    with cm as vcz1:
        append_function(vcz1, vcz2)


@click.command()
@click.argument("vcz", type=click.Path())
@click.argument("sample_id", type=str)
@impl
@zarr_backend_storage
def remove(vcz, sample_id, impl, zarr_backend_storage):
    """Remove a sample from vcz and overwrite with missing data"""
    if impl is None or impl == "zarr":
        from vczstore.zarr_impl import remove as remove_function
    elif impl == "cubed":
        from vczstore.cubed_impl import remove as remove_function
    elif impl == "xarray":
        from vczstore.xarray_impl import remove as remove_function
    else:
        raise ValueError(f"Unrecognised impl: {impl}")

    if zarr_backend_storage == "icechunk":
        from vczstore.icechunk_utils import icechunk_transaction

        cm = icechunk_transaction(vcz, "main", message="remove")
    else:
        cm = nullcontext(vcz)
    with cm as vcz:
        remove_function(vcz, sample_id)


@click.command()
@click.argument("vcz1", type=click.Path())
@click.argument("vcz2", type=click.Path())
def copy_store_to_icechunk(vcz1, vcz2):
    """Copy a Zarr store to a new Icechunk store"""
    from vczstore.icechunk_utils import (
        copy_store_to_icechunk as copy_store_to_icechunk_function,
    )

    copy_store_to_icechunk_function(vcz1, vcz2)


@click.group(cls=NaturalOrderGroup, name="vczstore")
def vczstore_main():
    pass


vczstore_main.add_command(append)
vczstore_main.add_command(remove)
vczstore_main.add_command(copy_store_to_icechunk)
