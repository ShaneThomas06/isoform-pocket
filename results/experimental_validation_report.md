# Experimental cavity validation

Each experimental PDB was reduced to its longest protein chain and all
HETATM records were removed before cavity detection. Experimental residue
numbers were sequence-mapped back to the corresponding AlphaFold model.
This validates monomeric cavities only, not oligomer-interface sites.

## RAC1/RAC1B splice-neighbourhood matches

| Isoform | AF cavity | AF volume (A3) | PDB | Experimental cavity | Experimental volume (A3) | Residue Jaccard | Inserted residues in AF cavity |
|---|---|---:|---|---|---:|---:|---:|
| RAC1 | KAB | 34.13 | 8S1N | KAE | 6.48 | 0.333 | 0 |

## Decision rule

A RAC1B cavity is retained only if it has meaningful experimental
residue overlap and does not depend primarily on the unresolved insertion.
The result still requires an independent pocket detector before being
described as a supported isoform-selective candidate.
