Initial experiments to try different technologies for vczlib to support the OFH Zarr POC.

### Operations
* Append to a Zarr store (atomically)
* Delete values from a Zarr store (atomically)
* Tag a Zarr store

### Matrix testing

* Transactions: none
* Distributed: single machine
* Zarr: v2

```shell
conda deactivate
conda env remove -n vczlib-poc-zarr-v2
conda create --name vczlib-poc-zarr-v2 -y 'python==3.12'
conda activate vczlib-poc-zarr-v2
pip install -e '.[dev]'
pip install -U -e ../vcztools  # delete-mask branch
pytest -vs
```

* Transactions: none
* Distributed: single machine
* Zarr: v3

```shell
conda deactivate
conda env remove -n vczlib-poc-zarr-v3
conda create --name vczlib-poc-zarr-v3 -y 'python==3.12'
conda activate vczlib-poc-zarr-v3
pip install -e '.[dev]'
pip install -U -e ../bio2zarr  # zarr-python-v3-fix branch
pip install -U -e ../vcztools  # delete-mask-zarr-python-v3-fix branch
pip install 'zarr>3'
pytest -vs
```