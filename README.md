Initial experiments to try different technologies for vczlib to support the OFH Zarr POC.

### Operations
* Append to a Zarr store (atomically)
* Delete values from a Zarr store (atomically)
* Tag a Zarr store

### Matrix testing

```shell
conda deactivate
conda env remove -n vczlib-poc-zarr-v2
conda create --name vczlib-poc-zarr-v2 -y 'python==3.12'
conda activate vczlib-poc-zarr-v2
pip install -e '.[dev]'
pip install -U -e ../vcztools  # delete-mask branch
pytest -vs -k "not icechunk"
```