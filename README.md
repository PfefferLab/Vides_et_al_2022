# Vides et. al 2022
Companion code for Vides et. al 2022.

This repository contains the code used to analyze the TIRF microscopy images in the Rab10-dependent LRRK2 recruitment in vitro assay.

The following steps were performed:
- LRRK2 recruitment movies were processed with the single molecule tracking and analysis framework TrackIT[[1]](#1). 
- Tracks of fluorescent intensity over time of individual LRRK2 spots were exported to csv files using a custom matlab script [tracks_export.m](tracks_export.m)
- Tracks were analyzed in R as reported in the [rab10_tirf.Rmd](Rmd/rab10_tirf.md) notebook.

## References
<a id="1">[1]</a> 
Kuhn, T., Hettich, J., Davtyan, R. et al. 
Single molecule tracking and analysis framework including theory-predicted parameter settings. Sci Rep 11, 9465 (2021). 
https://doi.org/10.1038/s41598-021-88802-7

