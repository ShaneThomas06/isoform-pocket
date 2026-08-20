# Seed-202 RAC1B/RAC1-p120 ARM-domain comparison

## Status

This report records the initial seed-202 observation and is retained for provenance. It is superseded by [`af3_multiseed_report.md`](af3_multiseed_report.md), which analyses seeds 202, 303, and 404 with five diffusion models per job.

## Seed-202 observation

The focused ARM-domain models had low global complex confidence. Mean ipTM was 0.168 for RAC1B and 0.190 for RAC1. Within this seed, RAC1B insertion residues 76-80 contacted CTNND1 positions 395, 430-437, and 474-477 more strongly than the matched RAC1 model.

| Construct | Mean ipTM | Mean pTM | Pairs >= 0.30 | Pairs >= 0.50 in all models |
|---|---:|---:|---:|---:|
| RAC1B | 0.168 | 0.640 | 15 | 0 |
| RAC1 | 0.190 | 0.652 | 6 | 0 |

## Multi-seed interpretation

The apparent RAC1B specificity did not reproduce at seeds 303 and 404. All six seed-level mean ipTM values were below 0.23, and no stable isoform-exclusive complex was recovered. The selected CTNND1 surface overlaps the established E-cadherin-binding groove in experimental structure 3L6Y and is at least 58.1 Angstrom from the implicated CTNND1 622-628 motif.

The seed-202 contact pattern must not be interpreted as a RAC1B-specific interface. See the multi-seed report and p120 structural-mapping report for the final analysis.
