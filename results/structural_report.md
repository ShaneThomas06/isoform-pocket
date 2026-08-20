# Structural comparison report

Structures were aligned only on sequence-identical residues using a Kabsch
least-squares superposition. The splice-altered residues therefore did not
determine the alignment frame.

## Global comparisons

| Gene | Reference | Compared model | Type | Mapped residues | RMSD (A) | P90 displacement (A) |
|---|---|---|---|---:|---:|---:|
| PKM | AF:PKM2 | AF:PKM1 | alphafold_isoform_pair | 486 | 0.228 | 0.283 |
| RAC1 | AF:RAC1 | AF:RAC1B | alphafold_isoform_pair | 192 | 0.982 | 0.535 |
| PKM | PDB:3SRF | AF:PKM1 | experimental_validation | 519 | 1.023 | 1.088 |
| PKM | PDB:3U2Z | AF:PKM2 | experimental_validation | 519 | 0.912 | 1.361 |
| RAC1 | PDB:8S1N | AF:RAC1 | experimental_validation | 172 | 1.152 | 0.987 |
| RAC1 | PDB:1RYF | AF:RAC1B | experimental_validation | 168 | 1.296 | 1.021 |
| RAC1 | PDB:1RYH | AF:RAC1B | experimental_validation | 169 | 1.351 | 0.94 |

## Isoform-local displacement after global superposition

- PKM reference residues 389-433: mean 0.378 A; maximum 0.968 A; n=45 identical residue pairs.
- RAC1 residues 66-75, immediately upstream of the insertion: mean 1.939 A; maximum 10.364 A; n=10.
- RAC1 residues 76-85, mapped to RAC1B residues 95-104 downstream of the insertion: mean 0.29 A; maximum 0.517 A; n=10.

## Interpretation boundary

A displacement is a geometric observation, not proof of a druggable site.
The next stage must test whether any displacement changes a reproducible
cavity or known interaction surface and whether experimental structures
support the same conclusion.
