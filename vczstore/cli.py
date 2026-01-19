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


@click.command()
@click.argument("vcz1", type=click.Path())
@click.argument("vcz2", type=click.Path())
@impl
def append(vcz1, vcz2, impl):
    """Append vcz2 to vcz1 in place"""
    if impl is None or impl == "zarr":
        from vczstore.zarr_impl import append as append_function
    elif impl == "cubed":
        from vczstore.cubed_impl import append as append_function
    elif impl == "xarray":
        from vczstore.xarray_impl import append as append_function
    else:
        raise ValueError(f"Unrecognised impl: {impl}")
    append_function(vcz1, vcz2)


@click.command()
@click.argument("vcz", type=click.Path())
@click.argument("sample_id", type=str)
@impl
def remove(vcz, sample_id, impl):
    """Remove a sample from vcz and overwrite with missing data"""
    if impl is None or impl == "zarr":
        from vczstore.zarr_impl import remove as remove_function
    elif impl == "cubed":
        from vczstore.cubed_impl import remove as remove_function
    elif impl == "xarray":
        from vczstore.xarray_impl import remove as remove_function
    else:
        raise ValueError(f"Unrecognised impl: {impl}")
    remove_function(vcz, sample_id)


@click.group(cls=NaturalOrderGroup, name="vczstore")
def vczstore_main():
    pass


vczstore_main.add_command(append)
vczstore_main.add_command(remove)
