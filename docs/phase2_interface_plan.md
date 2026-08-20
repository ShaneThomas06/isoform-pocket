# Phase 2: RAC1B interaction-surface analysis

## Decision

The primary partner will be **p120-catenin (CTNND1)**. **SmgGDS
(RAP1GDS1)** will be retained as a comparison partner if the first analysis is
successful.

## Why p120-catenin is first

The selection is based on direct experimental evidence, not simply a database
interaction:

- RAC1B binds p120-catenin more effectively than RAC1.
- The isolated 19-residue RAC1B insertion binds p120-catenin.
- Full-length RAC1B can compete with the isolated insertion for binding.
- The proposed recognition regions are flexible: RAC1B residues around 70-90
  and a p120-catenin region containing a polylysine motif around residues
  607-644.

This makes p120-catenin suitable for testing a disorder-to-order recognition
hypothesis. It does **not** establish one unique bound conformation.

Primary evidence:
[Orlichenko et al., Journal of Biological Chemistry (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2885194/)

Reviewed sequence records:
[CTNND1/O60716](https://www.uniprot.org/uniprotkb/O60716/entry) and
[RAP1GDS1/P52306](https://www.uniprot.org/uniprotkb/P52306/entry)

## Focused question

Does the RAC1B insertion create an isoform-specific physicochemical surface
that can recognize the flexible p120-catenin segment, even though it does not
form a stable small-molecule pocket by itself?

## Pre-registered hypotheses

### H1: insertion-driven recognition

The RAC1B 19-residue insertion contributes a surface chemistry pattern absent
from RAC1 and is the main source of isoform selectivity.

### H2: neighbouring-switch contribution

Residues 70-75 and the adjacent switch-II region contribute alongside the
insertion because their conformation or accessibility changes in RAC1B.

### H3: no single stable complex

Complex modelling will produce multiple plausible poses or locally low
confidence rather than one high-confidence rigid interface. This would be
consistent with flexible, disorder-to-order binding and must not be described
as proof of a binding mode.

## Minimal analysis

1. Retrieve and checksum the relevant CTNND1 and RAP1GDS1 isoform sequences.
2. Calculate charge, hydrophobicity, disorder propensity, and residue-level
   solvent exposure for the RAC1/RAC1B splice neighbourhood.
3. Define the exact p120-catenin construct used in the experimental study so
   its numbering is not confused with the many CTNND1 isoforms.
4. Generate several complex hypotheses, retaining all independent runs rather
   than selecting one attractive pose.
5. Compare interface-residue recurrence and confidence across runs.
6. Rank residue substitutions for experimental testing.

## Guardrails

- Do not force the low-confidence RAC1B insertion into one conformation.
- Do not call an AlphaFold complex prediction an experimentally validated
  interaction.
- Do not present transient RAC1B pockets as previously unreported: a 2025 study already
  reported molecular-dynamics-derived cryptic pockets.
- Report inconsistent complex poses as uncertainty, not as a result to hide.

## Intended deliverables

- an isoform-aware interface feature table;
- a multi-run interface recurrence plot;
- three to five candidate RAC1B residues for mutagenesis; and
- a concise hypothesis explaining why p120-catenin recognizes RAC1B more
  strongly than RAC1.

## Current status

The sequence-feature milestone is complete. Canonical CTNND1 numbering was
retained for p120 isoform 4, and the first mutation panel is:

1. D83A;
2. E78A;
3. D90A;
4. Y80A; and
5. K82A.

The acidic trio directly tests recognition of p120's polybasic loop. Y80A tests
an aromatic contribution, and K82A begins testing the insertion's basic
cluster. These are experimental hypotheses, not predicted contact residues.

## Multi-seed and experimental-structure result

The multi-run analysis is complete for matched RAC1 and RAC1B complexes at
seeds 202, 303, and 404, retaining five diffusion models per job.

- Every global pose remains low-confidence (mean ipTM below 0.23).
- The RAC1B-specific seed-202 pattern did not reproduce.
- A shared CTNND1 395-477 surface recurred, but its signal varied by seed.
- The surface geometry agrees with crystal structure 3L6Y (R^2 = 0.998).
- It lies at least 58.1 A from the experimentally implicated 622-628 motif.
- Six of seven selected residues contact E-cadherin in 3L6Y.

The evidence therefore does not support a RAC1B-specific binding site. Instead,
the study identifies a failure mode in which AlphaFold complex modelling can
place a partner against an existing structured peptide-binding groove when the
experimentally required recognition loop is disordered. This supports H3 and
the preregistered null-compatible guardrails.
