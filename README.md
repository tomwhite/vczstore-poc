# vczstore

Initial experiments to try different technologies for a vczstore POC.

### Operations

The key operations we want to support are:

* Append to a Zarr store (atomically)
* Remove values from a Zarr store (atomically)
* Tag a Zarr store

### Support matrix

Using zarr-python directly (single machine)

| zarr-python (format) | v2 (v2)            | v3 (v2)            | v3 (v3)            |
|----------------------|--------------------|--------------------|--------------------|
| **no transactions**  | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| **icechunk**         | N/A                | N/A                | :white_check_mark: |

Using Cubed (distributed)

| zarr-python (format) | v2 (v2)            | v3 (v2) | v3 (v3)            |
|----------------------|--------------------|---------|--------------------|
| **no transactions**  | :white_check_mark: | :x:     | :x:                |
| **icechunk**         | N/A                | N/A     | :white_check_mark: |

In addition the following features are not yet supported:

* Error checking (e.g. that the store and the vcz being appended have compatible fields)
* Allele harmonisation
* Cloud stores

### Demo

* Transactions: none
* Distributed: single machine
* Zarr: v2

```shell
conda activate vczstore-poc-zarr-v2

# Create some VCZ data
rm -rf data
mkdir data
vcf2zarr convert tests/data/vcf/sample-part1.vcf.gz data/store.vcz
vcf2zarr convert tests/data/vcf/sample-part2.vcf.gz data/sample-part2.vcf.vcz

# Show the samples in each
vcztools query -l data/store.vcz
vcztools query -l data/sample-part2.vcf.vcz

# Append data to the store
vczstore append data/store.vcz data/sample-part2.vcf.vcz
vcztools query -l data/store.vcz

# Remove a sample from the store
vczstore remove data/store.vcz NA00002
vcztools query -l data/store.vcz
```

* Transactions: none
* Distributed: cubed
* Zarr: v2

```shell
conda activate vczstore-poc-cubed-zarr-v2

# Create some VCZ data
rm -rf data
mkdir data
vcf2zarr convert tests/data/vcf/sample-part1.vcf.gz data/store.vcz
vcf2zarr convert tests/data/vcf/sample-part2.vcf.gz data/sample-part2.vcf.vcz

# Show the samples in each
vcztools query -l data/store.vcz
vcztools query -l data/sample-part2.vcf.vcz

# Append data to the store
vczstore append --impl cubed data/store.vcz data/sample-part2.vcf.vcz
vcztools query -l data/store.vcz

# Remove a sample from the store
vczstore remove --impl cubed data/store.vcz NA00002
vcztools query -l data/store.vcz
```

* Transactions: icechunk
* Distributed: single machine
* Zarr: v3, format 3

```shell
conda activate vczstore-poc-icechunk

# Create some VCZ data
rm -rf data
mkdir data
BIO2ZARR_ZARR_FORMAT=3 vcf2zarr convert tests/data/vcf/sample-part1.vcf.gz data/sample-part1.vcf.vcz
BIO2ZARR_ZARR_FORMAT=3 vcf2zarr convert tests/data/vcf/sample-part2.vcf.gz data/sample-part2.vcf.vcz

# Copy first vcz to an icechunk store
vczstore copy-store-to-icechunk data/sample-part1.vcf.vcz data/store.vcz
rm -rf data/sample-part1.vcf.vcz

# Show the samples in each
vcztools query -l data/store.vcz --zarr-backend-storage icechunk
vcztools query -l data/sample-part2.vcf.vcz

# Append data to the store
vczstore append data/store.vcz data/sample-part2.vcf.vcz --zarr-backend-storage icechunk
vcztools query -l data/store.vcz --zarr-backend-storage icechunk

# Remove a sample from the store
vczstore remove data/store.vcz NA00002 --zarr-backend-storage icechunk
vcztools query -l data/store.vcz --zarr-backend-storage icechunk
```

### Matrix testing

All (quick)

```shell
conda activate vczstore-poc-zarr-v2
pytest -vs

conda activate vczstore-poc-zarr-v3
pytest -vs

conda activate vczstore-poc-zarr-v3-f3
BIO2ZARR_ZARR_FORMAT=3 pytest -vs

conda activate vczstore-poc-cubed-zarr-v2
pytest -vs -k cubed

conda activate vczstore-poc-xarray-zarr-v2
pytest -vs -k xarray

conda activate vczstore-poc-icechunk
BIO2ZARR_ZARR_FORMAT=3 pytest -vs -k icechunk

conda activate vczstore-poc-icechunk-cubed
BIO2ZARR_ZARR_FORMAT=3 pytest -vs -k icechunk_cubed
```

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

* Transactions: icechunk
* Distributed: cubed
* Zarr: v3, format 3

```shell
conda deactivate
conda env remove -n vczstore-poc-icechunk-cubed
conda create --name vczstore-poc-icechunk-cubed -y 'python==3.12'
conda activate vczstore-poc-icechunk-cubed
pip install -e '.[dev]'
pip install -U 'git+https://github.com/sgkit-dev/bio2zarr.git'
pip install -U 'git+https://github.com/tomwhite/vcztools.git@sample-mask-icechunk'
# pip install -U -e ../bio2zarr  # zarr-format-3 branch
# pip install -U -e ../vcztools  # sample-mask-icechunk branch
pip install 'zarr>3' icechunk hypothesis
pip install -U 'git+https://github.com/cubed-dev/cubed.git@update'
# pip install -e '../cubed[diagnostics]' # update branch
pip install rich
BIO2ZARR_ZARR_FORMAT=3 pytest -vs -k icechunk_cubed
```