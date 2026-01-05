Initial experiments to try different technologies for vczstore to support the OFH Zarr POC.

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
conda env remove -n vczstore-poc-zarr-v2
conda create --name vczstore-poc-zarr-v2 -y 'python==3.12'
conda activate vczstore-poc-zarr-v2
pip install -e '.[dev]'
pip install -U 'git+https://github.com/tomwhite/vcztools.git@sample-mask'
# pip install -U -e ../vcztools  # sample-mask branch
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

* Transactions: none
* Distributed: single machine
* Zarr: v3, format 3

```shell
conda deactivate
conda env remove -n vczlib-poc-zarr-v3-f3
conda create --name vczlib-poc-zarr-v3-f3 -y 'python==3.12'
conda activate vczlib-poc-zarr-v3-f3
pip install -e '.[dev]'
pip install -U -e ../bio2zarr  # zarr-python-v3-fix-zarr-format-3 branch
pip install -U -e ../vcztools  # delete-mask-zarr-python-v3-fix branch
pip install 'zarr>3'
pytest -vs
```

* Transactions: none
* Distributed: cubed
* Zarr: v2

```shell
conda deactivate
conda env remove -n vczlib-poc-cubed-zarr-v2
conda create --name vczlib-poc-cubed-zarr-v2 -y 'python==3.12'
conda activate vczlib-poc-cubed-zarr-v2
pip install -e '.[dev]'
pip install -U -e ../vcztools  # delete-mask branch
pip install -e '../cubed[diagnostics]' # update branch
pytest -vs -k cubed
```

* Transactions: icechunk
* Distributed: single machine
* Zarr: v3, format 3

```shell
conda deactivate
conda env remove -n vczlib-poc-icechunk
conda create --name vczlib-poc-icechunk -y 'python==3.12'
conda activate vczlib-poc-icechunk
pip install -e '.[dev]'
pip install -U -e ../bio2zarr  # zarr-python-v3-fix-zarr-format-3 branch
pip install -U -e ../vcztools  # delete-mask-zarr-python-v3-fix-icechunk branch
pip install 'zarr>3' icechunk
pytest -vs -k icechunk
```