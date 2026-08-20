# Milestone 1: Is there a fixed RAC1B-specific pocket?

## Short answer

No reproducible fixed pocket was found beside the RAC1B insertion.

This is not a failed analysis. It is an evidence-based negative result that
prevents the project from making an attractive but unsupported druggability
claim.

## What was compared

The project used PKM1/PKM2 as a positive control and RAC1/RAC1B as the discovery
pair. It combined:

1. reviewed UniProt isoform sequences;
2. AlphaFold models and per-residue confidence;
3. experimental crystal structures;
4. pyKVFinder geometric cavity detection; and
5. P2Rank machine-learning pocket prediction.

## Key observations

### PKM positive control

PKM1 and PKM2 have different residues at positions 389-433, but both
isoform-specific regions have very high AlphaFold confidence. The two predicted
cores superimpose closely (0.228 Å RMSD), so the workflow correctly recovers a
localized, confidently modelled isoform change rather than a different global
fold.

### RAC1B insertion

RAC1B contains a 19-residue insertion at positions 76-94. Its mean AlphaFold
confidence is only 45.27, and only 5.3% of inserted residues reach pLDDT 70.
The insertion should therefore be interpreted as flexible or structurally
uncertain. not as one stable conformation.

The adjacent RAC1/RAC1B cores remain similar (0.982 Å RMSD). Most of the
apparent upstream displacement is concentrated at residues 74-75 next to the
insertion, while the downstream segment aligns closely.

## Pocket result

pyKVFinder detected a large 355.54 Å³ cavity involving 16 low-confidence
inserted residues. P2Rank did not reproduce a splice-neighbourhood pocket, and
neither RAC1B crystal structure supported it. The cavity is therefore most
consistent with a model-conformation artefact.

pyKVFinder also matched a smaller cavity between RAC1 and RAC1B AlphaFold
models, but this cavity again lacked support in either RAC1B crystal structure.
P2Rank's highest-ranked RAC1B pocket was not located at the splice
neighbourhood.

## Defensible conclusion

The available evidence does not support docking compounds into a fixed,
RAC1B-selective pocket beside the insertion. The stronger biological hypothesis
is that the flexible insertion changes transient interaction surfaces,
conformational dynamics, or partner recognition.

## Recommended phase 2

The next phase will test whether RAC1B changes protein-partner recognition:

1. define a focused RAC1/RAC1B surface neighbourhood around the insertion;
2. map experimentally reported interaction residues where available;
3. compare solvent exposure, charge, and hydrophobicity between isoforms;
4. model selected protein complexes cautiously, using confidence and
   experimental evidence as filters; and
5. nominate testable interface mutations rather than speculative ligands.

The deliverable will be a ranked table of interface residues and a small set of
experimentally testable hypotheses.
