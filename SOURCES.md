# Sources, mechanisms and limitations

| Modeled mechanism | Source | Status in this simulator |
|---|---|---|
| Power depends on donors, cells per donor and depth with diminishing returns | Schmid KT et al. scPower accelerates and optimizes the design of multi-sample single cell transcriptomic studies. Nat Commun 12:6625 (2021) | Sourced; saturation curve 1 - exp(-reads/8000) is an invented fixture |
| Batch and biological effects are not separable when conditions are confounded with batch | Song F, Chan GMA, Wei Y. Flexible experimental designs for valid single-cell RNA-sequencing experiments allowing batch effects correction. Nat Commun 11:3274 (2020) | Sourced design principle |
| Treating cells as replicates inflates false discoveries; donor-level (pseudobulk) analysis controls them | Squair JW et al. Confronting false discoveries in single-cell differential expression. Nat Commun 12:5692 (2021) | Sourced; reproduced by the `cell_level` comparison route |
| UMI counts are well described by a negative binomial without zero inflation | Svensson V. Droplet scRNA-seq is not zero-inflated. Nat Biotechnol 38:147-150 (2020) | Sourced; dispersion 0.15 is invented |
| Ambient RNA contaminates some libraries | Young MD, Behjati S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. GigaScience 9:giaa151 (2020) | Sourced qualitatively; contamination fractions are invented |

Invented fixtures: the 40-gene panel, marker levels, donor and batch effect sizes, effect magnitudes and costs.

Limitations: a small panel instead of a transcriptome, perfect marker separation of three types, no doublets, no cell-cycle or continuous states, independent donors, and costs that are illustrative only.
