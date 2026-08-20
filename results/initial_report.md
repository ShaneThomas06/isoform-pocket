# Initial input and confidence report

## What was checked

All four reviewed isoform sequences and AlphaFold models were retrieved,
checksummed, and compared. Five experimental structures were also retrieved
for the next structural-validation stage.

## AlphaFold model quality

| Protein | Length | Mean pLDDT | Fraction >=70 | Entry |
|---|---:|---:|---:|---|
| PKM2 | 531 | 96.84 | 98.7% | AF-P14618-F1 |
| PKM1 | 531 | 96.03 | 98.5% | AF-P14618-2-F1 |
| RAC1 | 192 | 93.84 | 93.8% | AF-P63000-F1 |
| RAC1B | 211 | 88.14 | 83.9% | AF-P63000-2-F1 |

## Splice differences

- **PKM2 vs PKM1:** replacement; reference 389-433, comparator 389-433.
- **RAC1 vs RAC1B:** insert; reference 76-75, comparator 76-94.

## Region-level confidence

| Protein | Region | Residues | Mean pLDDT | Minimum | Fraction >=70 |
|---|---|---:|---:|---:|---:|
| PKM2 | alternative_exon_region | 389-433 | 98.37 | 96.94 | 100.0% |
| PKM1 | alternative_exon_region | 389-433 | 94.63 | 83.75 | 100.0% |
| RAC1 | insertion_site_neighborhood | 66-85 | 97.88 | 94.69 | 100.0% |
| RAC1B | upstream_flank | 66-75 | 78.53 | 43.38 | 70.0% |
| RAC1B | inserted_segment | 76-94 | 45.27 | 35.34 | 5.3% |
| RAC1B | downstream_flank | 95-104 | 97.7 | 91.19 | 100.0% |

## Decision

PKM1/PKM2 passes the fixed-structure control gate. The RAC1B inserted
segment is treated as flexible unless experimental evidence supports a
specific bound conformation. Confident neighbouring regions proceed to
structural and interface comparison.
