# Consensus pocket decision

## P2Rank results for RAC1 and RAC1B

| Isoform | Pocket | Rank | Score | Probability | Residues | Touches splice neighbourhood |
|---|---|---:|---:|---:|---|---|
| RAC1 | pocket1 | 1 | 6.35 | 0.266 | 13;15;17;18;28;32;115;116;118;119;159 | False |
| RAC1 | pocket2 | 2 | 5.22 | 0.2 | 12;13;15;16;17;18;31;32;33;34;35;60;61 | False |
| RAC1B | pocket1 | 1 | 10.8 | 0.498 | 12;13;14;15;16;17;18;19;28;29;32;33;34;35;60;61;134;135;137;138;178;179 | False |

## Cross-method decision

- pyKVFinder detected 3 RAC1/RAC1B cavities touching the splice neighbourhood.
- P2Rank detected no pocket touching the RAC1B insertion or immediate flanks.
- The largest pyKVFinder RAC1B cavity depends on 16 low-confidence inserted residues.
- The smaller high-confidence model cavity was not recovered in either RAC1B crystal structure.

## Current conclusion

There is currently **no reproducible evidence for an isoform-selective, fixed
small-molecule pocket near the RAC1B insertion**. The results instead
support treating this region as a dynamic interaction surface. This is a
scientifically useful negative result and prevents a false druggability claim.

The next discovery analysis should test altered partner-binding surfaces
and disorder-to-order behaviour rather than perform docking into the
unsupported AlphaFold cavity.
