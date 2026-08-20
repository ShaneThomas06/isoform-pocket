# Study design

## Primary question

Does the RAC1B-specific insertion alter an adjacent, structurally supported
pocket or protein-interaction surface relative to canonical RAC1?

## Systems

### Positive control: PKM2 versus PKM1

PKM2 (`P14618-1`, canonical AlphaFold query `P14618`) and PKM1 (`P14618-2`)
have the same length but replace a 45-residue region. The control tests whether
the workflow correctly localizes a known isoform-dependent change.

Experimental controls:

- `3SRF`: human PKM1.
- `3U2Z`: human PKM2 with an allosteric activator.

### Discovery system: RAC1 versus RAC1B

RAC1B (`P63000-2`) contains a 19-residue insertion after RAC1 residue 75. The
insertion is flexible, so fixed-pocket analysis is restricted to confident
neighbouring residues. Interface and disorder-to-order hypotheses remain in
scope.

Experimental controls:

- `8S1N`: human RAC1-GDP.
- `1RYF`: human RAC1B-GDP.
- `1RYH`: human RAC1B-GppNHp.

## Preregistered hypotheses

- **H1:** PKM1 and PKM2 retain a conserved global fold while their strongest
  localized sequence/structure differences occur around residues 389-433.
- **H2:** RAC1B changes the conformation or accessibility of at least one
  confident residue group adjacent to its insertion.
- **H0:** No reproducible RAC1B-specific pocket or interface remains after confidence
  filtering and experimental comparison. This is an acceptable result.

## Evidence gates

A structural feature is reportable only when:

1. sequence and residue mappings are unambiguous;
2. surrounding local structure has mean pLDDT >= 70, or experimental support;
3. it is reproducible across reasonable analysis settings;
4. its interpretation does not depend on one arbitrary low-confidence loop;
5. limitations are reported alongside its score.

Strong pocket claims will eventually require agreement between two independent
pocket detectors.

## Claim boundaries

The study may identify or prioritize a candidate pocket/interface. It will not
claim experimental binding, measured affinity, clinical significance, or a
validated drug target.
