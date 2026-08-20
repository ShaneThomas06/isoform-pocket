# RAC1B-p120 interface feature report

## Numbering and evidence

The binding experiments used short p120-catenin isoform 4. This analysis
retains canonical human CTNND1/O60716 coordinates for the shared ARM-domain
sequence (607-644). RAC1B uses P63000-2 coordinates. No complex structure
is assumed.

## Sequence-level comparison

| Region | Sequence | Net charge | Positive | Negative | Hydrophobic | Mean hydropathy |
|---|---|---:|---:|---:|---:|---:|
| RAC1B 76-94 | `VGETYGKDITSRGKDKPIA` | 1 | 21.1% | 15.8% | 26.3% | -0.947 |
| p120 607-644 | `AERYQEAAPNVANNTGPHAASCFGAKKGKDEWFSRGKK` | 3 | 18.4% | 10.5% | 31.6% | -1.161 |

Both segments have mixed charge and flexible residues. Because the p120
segment is net positive, RAC1B acidic residues are direct experimental
candidates, but composition alone cannot specify a binding pose.

## First-pass mutagenesis candidates

| Rank | RAC1B residue | Suggested test | Reason |
|---:|---|---|---|
| 1 | D83 | D83A | tests recognition of the p120 polybasic loop |
| 2 | E78 | E78A | tests recognition of the p120 polybasic loop |
| 3 | D90 | D90A | tests recognition of the p120 polybasic loop |
| 4 | Y80 | Y80A | tests a possible transient aromatic anchor |
| 5 | K82 | K82A | tests the insertion's basic cluster |

## Interpretation

These candidates test electrostatic, aromatic, and conformational
features of the insertion. The ranking is not a prediction of
experimental effect. Because insertion pLDDT is low, AlphaFold exposure
values are exploratory metadata and do not define a fixed interface.
