# bcftools view tests/data/vcf/alleles-1.vcf -O z > tests/data/vcf/alleles-1.vcf.gz
# bcftools view tests/data/vcf/alleles-2.vcf -O z > tests/data/vcf/alleles-2.vcf.gz
# bcftools index -c tests/data/vcf/alleles-1.vcf.gz
# bcftools index -c tests/data/vcf/alleles-2.vcf.gz

# merging
# bcftools merge tests/data/vcf/alleles-1.vcf.gz tests/data/vcf/alleles-2.vcf.gz
# bcftools merge tests/data/vcf/alleles-1.vcf.gz tests/data/vcf/alleles-2.vcf.gz -O z \
#  > tests/data/vcf/alleles-merged.vcf.gz

import numpy as np
import zarr

from vczstore.allele_harmonisation import harmonise_alleles
from vczstore.zarr_impl import append_harmonise

from .utils import (
    compare_vcf_and_vcz,
    convert_vcf_to_vcz,
)


def test_append_different_alleles(tmp_path):
    print(tmp_path)

    vcz1 = convert_vcf_to_vcz("alleles-1.vcf.gz", tmp_path)
    vcz2 = convert_vcf_to_vcz("alleles-2.vcf.gz", tmp_path)

    root1 = zarr.open(vcz1, mode="r")
    root2 = zarr.open(vcz2, mode="r")
    print("variant_allele 1", root1["variant_allele"][:])
    print("variant_allele 2", root2["variant_allele"][:])

    variant_allele1 = root1["variant_allele"][:]
    variant_allele2 = root2["variant_allele"][:]

    variant_allele1_updated, variant_allele2_mapping = harmonise_alleles(
        variant_allele1, variant_allele2
    )

    np.testing.assert_array_equal(variant_allele1_updated, [["A", "C", "G"]])
    np.testing.assert_array_equal(variant_allele2_mapping, [[0, 2]])


def test_append_harmonise(tmp_path):
    print(tmp_path)

    vcz1 = convert_vcf_to_vcz("alleles-1.vcf.gz", tmp_path)
    vcz2 = convert_vcf_to_vcz("alleles-2.vcf.gz", tmp_path)

    append_harmonise(vcz1, vcz2)

    # check equivalence with original VCF
    compare_vcf_and_vcz(
        tmp_path,
        "view --no-version",
        "alleles-merged.vcf.gz",
        "view --no-version",
        vcz1,
    )
