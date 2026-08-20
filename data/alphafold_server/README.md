# AlphaFold Server input directory

AlphaFold Server result archives are not committed because the extracted jobs
occupy about 600 MB. Place each extracted result under the directory name shown
below before running the multi-seed analyses:

```text
RAC1_p120ARM350-824_seed202/
RAC1_p120ARM350-824_seed303/
RAC1_p120ARM350-824_seed404/
RAC1B_p120ARM350-824_seed202/
RAC1B_p120ARM350-824_seed303/
RAC1B_p120ARM350-824_seed404/
```

Each directory must contain an `extracted/` folder with the five
`full_data`, five `summary_confidences`, and five model CIF files supplied by
AlphaFold Server. The analysis scripts check that all 30 models are present.

The downloaded files remain subject to the AlphaFold Server terms included in
each result archive.
