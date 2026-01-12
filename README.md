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
conda env remove -n vczstore-poc-zarr-v3
conda create --name vczstore-poc-zarr-v3 -y 'python==3.12'
conda activate vczstore-poc-zarr-v3
pip install -e '.[dev]'
pip install -U 'git+https://github.com/sgkit-dev/bio2zarr.git'
pip install -U 'git+https://github.com/tomwhite/vcztools.git@sample-mask'
# pip install -U -e ../bio2zarr  # main branch
# pip install -U -e ../vcztools  # sample-mask branch
pip install 'zarr>3'
pytest -vs
```

* Transactions: none
* Distributed: single machine
* Zarr: v3, format 3

```shell
conda deactivate
conda env remove -n vczstore-poc-zarr-v3-f3
conda create --name vczstore-poc-zarr-v3-f3 -y 'python==3.12'
conda activate vczstore-poc-zarr-v3-f3
pip install -e '.[dev]'
pip install -U 'git+https://github.com/sgkit-dev/bio2zarr.git'
pip install -U 'git+https://github.com/tomwhite/vcztools.git@sample-mask'
# pip install -U -e ../bio2zarr  # zarr-format-3 branch
# pip install -U -e ../vcztools  # sample-mask branch
pip install 'zarr>3'
BIO2ZARR_ZARR_FORMAT=3 pytest -vs
```

* Transactions: none
* Distributed: cubed
* Zarr: v2

```shell
conda deactivate
conda env remove -n vczstore-poc-cubed-zarr-v2
conda create --name vczstore-poc-cubed-zarr-v2 -y 'python==3.12'
conda activate vczstore-poc-cubed-zarr-v2
pip install -e '.[dev]'
pip install -U 'git+https://github.com/tomwhite/vcztools.git@sample-mask'
# pip install -U -e ../vcztools  # sample-mask branch
pip install -U 'git+https://github.com/cubed-dev/cubed.git@update'
# pip install -e '../cubed[diagnostics]' # update branch
pip install rich
pytest -vs -k cubed
```

* Transactions: none
* Distributed: xarray
* Zarr: v2

```shell
conda deactivate
conda env remove -n vczstore-poc-xarray-zarr-v2
conda create --name vczstore-poc-xarray-zarr-v2 -y 'python==3.12'
conda activate vczstore-poc-xarray-zarr-v2
pip install -e '.[dev]' xarray
pip install -U 'git+https://github.com/tomwhite/vcztools.git@sample-mask'
# pip install -U -e ../vcztools  # sample-mask branch
pytest -vs -k xarray
```

* Transactions: icechunk
* Distributed: single machine
* Zarr: v3, format 3

```shell
conda deactivate
conda env remove -n vczstore-poc-icechunk
conda create --name vczstore-poc-icechunk -y 'python==3.12'
conda activate vczstore-poc-icechunk
pip install -e '.[dev]'
pip install -U 'git+https://github.com/sgkit-dev/bio2zarr.git'
pip install -U 'git+https://github.com/tomwhite/vcztools.git@sample-mask-icechunk'
# pip install -U -e ../bio2zarr  # zarr-format-3 branch
# pip install -U -e ../vcztools  # sample-mask-icechunk branch
pip install 'zarr>3' icechunk hypothesis
BIO2ZARR_ZARR_FORMAT=3 pytest -vs -k icechunk
```