Initial experiments to try different technologies for vczlib to support the OFH Zarr POC.

Operations
* Append to a Zarr store (atomically)
* Delete values from a Zarr store (atomically)
* Tag a Zarr store

This shows how we can remove a sample from a Zarr store (zarr-python/format v2 currently)

```shell
conda activate vczlib-poc-zarr-v2
rm -rf sample.vcf.vcz; cp -r ~/workspace/vcztools/vcz_test_cache/sample.vcf.vcz sample.vcf.vcz
pytest -vs tests/test_remove.py

# then in vcztools dir, delete-mask branch, vcztools-3.12 venv
python -m vcztools view ~/workspace/vczlib-poc/sample.vcf.vcz
```