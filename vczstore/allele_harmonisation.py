import numpy as np
from vcztools.utils import search


def harmonise_alleles(variant_allele1, variant_allele2):
    n_variants = variant_allele1.shape[0]
    # TODO: check variant_allele2 has the same number of variants

    # first find new max alleles
    max_alt_alleles = 0
    for i in range(n_variants):
        alt1 = variant_allele1[i][1:]
        alt2 = variant_allele2[i][1:]
        # for each alt allele in alt2
        new_alleles = np.setdiff1d(alt2, alt1)
        max_alt_alleles = max(len(alt1) + len(new_alleles), max_alt_alleles)
    max_alleles = max_alt_alleles + 1

    variant_allele1_updated = np.empty(
        (n_variants, max_alleles), dtype=variant_allele1.dtype
    )
    variant_allele2_mapping = np.zeros(
        (n_variants, variant_allele2.shape[1]), dtype=np.int64
    )

    for i in range(n_variants):
        ref1 = variant_allele1[i][0]
        ref2 = variant_allele2[i][0]

        if ref1 != ref2:
            raise ValueError("References must be the same when appending")

        alt1 = variant_allele1[i][1:]
        alt2 = variant_allele2[i][1:]

        if np.all(alt1 == alt2):
            # TODO: need to pad with fill
            variant_allele1_updated[i] = variant_allele1[i]

        else:
            # for each alt allele in alt2
            new_alleles = np.setdiff1d(alt2, alt1)
            updated = np.append(variant_allele1[i], new_alleles, axis=0)
            # TODO: need to pad with fill
            variant_allele1_updated[i] = updated

            mapping = search(updated, variant_allele2)
            variant_allele2_mapping[i] = mapping

    return variant_allele1_updated, variant_allele2_mapping


def remap_gt(gt, variant_allele_mapping):
    num_variants = gt.shape[0]
    num_samples = gt.shape[1]
    ploidy = gt.shape[2]
    for i in range(num_variants):
        for j in range(num_samples):
            for k in range(ploidy):
                val = gt[i, j, k]
                if val >= 0:
                    gt[i, j, k] = variant_allele_mapping[i, val]
    return gt
