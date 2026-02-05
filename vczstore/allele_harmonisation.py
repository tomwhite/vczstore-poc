import numpy as np
from vcztools.utils import search


def harmonise_alleles(variant_allele, variant_allele_new):
    """
    Harmonise alleles at each variant site.

    Returns two arrays:

    ``variant_allele_updated``: a ``variant_allele`` array such that each
    variant site contains the union of the alleles in ``variant_allele``
    and ``variant_allele_new``

    ``variant_allele_new_mapping``: an int array with the same shape as
    ``variant_allele_new`` where each value is an index into the allele in
    ``variant_allele_updated``. This is similar to the local alleles
    ``LAA`` field in VCF 4.5 with an extra leading column (which is always
    zero since REF is unchanged).
    """
    n_variants = variant_allele.shape[0]
    n_variants_new = variant_allele_new.shape[0]

    if n_variants != n_variants_new:
        raise ValueError(
            f"Different number of variants: store has {n_variants} variants, "
            f"but appended data has {n_variants_new}"
        )

    # first find new max alleles
    max_alt_alleles = 0
    for i in range(n_variants):
        alt = variant_allele[i][1:]
        alt_new = variant_allele_new[i][1:]
        # for each alt allele in alt_new
        new_alleles = np.setdiff1d(alt_new, alt)
        max_alt_alleles = max(len(alt) + len(new_alleles), max_alt_alleles)
    max_alleles = max_alt_alleles + 1

    variant_allele_updated = np.full(
        (n_variants, max_alleles), fill_value="", dtype=variant_allele.dtype
    )
    variant_allele_new_mapping = np.full(
        (n_variants, variant_allele_new.shape[1]), fill_value=-2, dtype=np.int64
    )

    for i in range(n_variants):
        ref = variant_allele[i][0]
        ref_new = variant_allele_new[i][0]

        if ref != ref_new:
            # TODO: would be more useful to have chr/pos/id in message
            raise ValueError(
                f"References don't match at index {i}: store is '{ref}', "
                f"but appended data is '{ref_new}'"
            )

        alt = variant_allele[i][1:]
        alt_new = variant_allele_new[i][1:]

        # for each alt allele in alt_new
        new_alleles = np.setdiff1d(alt_new, alt)
        updated = np.append(variant_allele[i], new_alleles, axis=0)
        variant_allele_updated[i][: updated.shape[0]] = updated

        # remove fill values at end of array
        variant_allele_new_no_fill = variant_allele_new[i][variant_allele_new[i] != ""]
        mapping = search(updated, variant_allele_new_no_fill)
        variant_allele_new_mapping[i][: variant_allele_new_no_fill.shape[0]] = mapping

    return variant_allele_updated, variant_allele_new_mapping


def remap_gt(gt, variant_allele_new_mapping):
    """Update a genotype array in-place using a mapping from ``harmonise_alleles``."""
    num_variants = gt.shape[0]
    num_samples = gt.shape[1]
    ploidy = gt.shape[2]
    for i in range(num_variants):
        for j in range(num_samples):
            for k in range(ploidy):
                val = gt[i, j, k]
                if val >= 0:
                    gt[i, j, k] = variant_allele_new_mapping[i, val]
